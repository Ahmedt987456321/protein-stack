"""Step 09 (v0.2) - identity-aware evidence fusion with an InterPro domain stream.

Streams, per test protein and GO term:
  seq(t)  max identity among sequence neighbours annotated with t  (= v0.1 arm A)
  str(t)  max TM-score among structural neighbours annotated with t
  dom(t)  max over the protein's InterPro entries of P(t | domain) estimated on
          the training side (min support 3, floor 0.02), or 1.0 for terms in the
          curated interpro2go mapping (propagated)

Fusion (arm C), noisy-OR with BOTH auxiliary streams gated by how strong the
protein's sequence evidence is (s_seq = best identity to the training set):

  w_str    = alpha  if s_seq < tau    else beta
  w_dom    = gamma  if s_seq < tau_d  else delta
  score(t) = 1 - (1 - seq(t)) * (1 - w_str*str(t)) * (1 - w_dom*dom(t))

(tau, alpha, beta, gamma, tau_d, delta) are fitted on an internal validation
split of the training set (never on test), maximising the macro-average Fmax
over branch x identity-bin cells - which explicitly punishes high-identity
dilution.

Also writes single-stream arms S (structure-only) and D (domain-only) so the
report can verify that fusion beats every single stream.

Run with --searches-only to just produce the validation-vs-subtrain searches.

Outputs:
  results/pred_armC.tsv, pred_armS.tsv, pred_armD.tsv
  results/fusion_params.json
"""
import json
import random
import shutil
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import (
    data_dir,
    link_or_copy,
    load_config,
    read_fasta,
    results_dir,
    run_tool,
    write_fasta,
)
from pis.eval import fmax
from pis.go import GoDag
from pis.streams import (
    BIN_LABELS,
    BRANCHES,
    bin_of,
    domain_stats,
    domain_stream,
    fuse,
    load_domains,
    load_interpro2go,
    parse_raw_seq,
    parse_raw_struct,
    split_by_branch,
    transfer_stream,
)

MIN_CELL = 10  # ignore tiny validation cells when fitting

GRID_TAU = [0.35, 0.5, 2.0]   # 2.0 = never gate (structure always at alpha)
GRID_ALPHA = [0.7, 1.0]
GRID_BETA = [0.0, 0.4]
GRID_GAMMA = [0.4, 0.8, 1.0]
GRID_TAU_D = [0.5, 0.8, 2.0]  # 2.0 = domain stream always at gamma
GRID_DELTA = [0.0, 0.5]


def load_split(path):
    train, test = [], []
    with open(path, encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, split, _b = line.rstrip("\n").split(",")
            (train if split == "train" else test).append(acc)
    return sorted(train), sorted(test)


def load_annotations(path):
    ann = defaultdict(set)
    with open(path, encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, term, _aspect = line.rstrip("\n").split("\t")
            ann[acc].add(term)
    return ann


def load_clean_hits(path, allowed_targets):
    """data/*_hits.tsv (query, target, sim) -> {query: [(target, sim)]}."""
    hits = defaultdict(list)
    with open(path, encoding="utf-8") as f:
        next(f)
        for line in f:
            q, t, sim = line.rstrip("\n").split("\t")
            if t in allowed_targets and q != t:
                hits[q].append((t, float(sim)))
    return hits


def macro_fmax(cells, fused_by_branch, truth_by_branch):
    scores = []
    for branch, proteins in cells:
        f, _ = fmax(proteins, fused_by_branch[branch], truth_by_branch[branch])
        scores.append(f)
    return sum(scores) / len(scores) if scores else 0.0


def ensure_val_searches(cfg, d, val, subtrain, seqs):
    """Create validation fasta/structure inputs and run both searches (cached)."""
    if not (d / "val.fasta").exists():
        write_fasta(d / "val.fasta", {a: seqs[a] for a in val})
        write_fasta(d / "subtrain.fasta", {a: seqs[a] for a in subtrain})
    for name, members in (("struct_val", val), ("struct_subtrain", subtrain)):
        out = d / name
        if not out.exists():
            out.mkdir()
            for acc in members:
                src = d / "structures" / (acc + ".pdb")
                if src.exists():
                    link_or_copy(src, out / (acc + ".pdb"))

    th = str(cfg["tools"]["threads"])
    if not (d / "val_seq_hits_raw.tsv").exists():
        run_tool(cfg, "mmseqs", [
            "easy-search", d / "val.fasta", d / "subtrain.fasta",
            d / "val_seq_hits_raw.tsv", d / "tmp_val_seq",
            "--format-output", "query,target,pident,evalue",
            "-e", "10", "-s", "7.5",
            "--max-seqs", str(cfg["search"]["max_seqs"]), "--threads", th,
        ])
        shutil.rmtree(d / "tmp_val_seq", ignore_errors=True)
    if not (d / "val_struct_hits_raw.tsv").exists():
        run_tool(cfg, "foldseek", [
            "easy-search", d / "struct_val", d / "struct_subtrain",
            d / "val_struct_hits_raw.tsv", d / "tmp_val_struct",
            "--format-output", "query,target,alntmscore,bits,evalue",
            "-e", str(cfg["search"]["evalue"]),
            "--max-seqs", str(cfg["search"]["max_seqs"]), "--threads", th,
        ])
        shutil.rmtree(d / "tmp_val_struct", ignore_errors=True)


def write_preds(path, preds_by_branch):
    n = 0
    merged = defaultdict(dict)
    for per_acc in preds_by_branch.values():
        for acc, scores in per_acc.items():
            merged[acc].update(scores)
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("accession\tgo_term\tscore\n")
        for acc in sorted(merged):
            for term, score in sorted(merged[acc].items()):
                f.write("{}\t{}\t{:.4f}\n".format(acc, term, score))
                n += 1
    return n


def main():
    cfg = load_config()
    d = data_dir(cfg)
    r = results_dir(cfg)
    edges = cfg["split"]["bin_edges"]
    transfer_e = cfg["search"]["evalue"]

    dag = GoDag(d / "go-basic.obo")
    train, test = load_split(d / "split.csv")
    raw_ann = load_annotations(d / "annotations.tsv")
    seqs = read_fasta(d / "final.fasta")

    rng = random.Random(cfg["seed"] + 1)
    shuffled = list(train)
    rng.shuffle(shuffled)
    n_val = round(0.2 * len(shuffled))
    val = sorted(shuffled[:n_val])
    subtrain = sorted(shuffled[n_val:])
    print("Fusion fit: {} subtrain / {} validation".format(len(subtrain), len(val)))

    ensure_val_searches(cfg, d, val, subtrain, seqs)
    if "--searches-only" in sys.argv:
        print("Validation searches ready (--searches-only).")
        return

    prop = {acc: dag.propagate(raw_ann[acc]) for acc in raw_ann}
    domains = load_domains(d / "domains.tsv")
    ip2go = load_interpro2go(d / "interpro2go.txt", dag)
    print("Domains: {} proteins with entries; interpro2go: {} mapped entries".format(
        len(domains), len(ip2go)))

    # ---------------- validation-side streams (subtrain knowledge only) ----
    sub_set = set(subtrain)
    val_s_seq, val_seq_transfer = parse_raw_seq(d / "val_seq_hits_raw.tsv", sub_set, transfer_e)
    val_struct = parse_raw_struct(d / "val_struct_hits_raw.tsv", sub_set)
    sub_terms = {a: prop[a] for a in subtrain if a in prop}
    val_seq_stream = split_by_branch(transfer_stream(val_seq_transfer, sub_terms), dag)
    val_str_stream = split_by_branch(transfer_stream(val_struct, sub_terms), dag)
    sub_probs = domain_stats(subtrain, prop, domains)
    val_dom_stream = split_by_branch(domain_stream(val, domains, sub_probs, ip2go), dag)

    val_truth = {b: {} for b in BRANCHES}
    for acc in val:
        for b in BRANCHES:
            sub = {t for t in prop.get(acc, ()) if dag.branch(t) == b}
            if sub:
                val_truth[b][acc] = sub
    val_bins = {acc: bin_of(val_s_seq.get(acc, 0.0), edges) for acc in val}

    cells = []
    for b in BRANCHES:
        for lbl in BIN_LABELS:
            members = sorted(a for a in val_truth[b] if val_bins[a] == lbl)
            if len(members) >= MIN_CELL:
                cells.append((b, members))
    print("Fitting on {} validation cells".format(len(cells)))

    # ---------------- grid search ------------------------------------------
    best = None
    for tau in GRID_TAU:
        betas = GRID_BETA if tau <= 1.0 else [0.0]      # beta unused when never gated
        for alpha in GRID_ALPHA:
            for beta in betas:
                for gamma in GRID_GAMMA:
                    for tau_d in GRID_TAU_D:
                        deltas = GRID_DELTA if tau_d <= 1.0 else [0.0]
                        for delta in deltas:
                            fused = {
                                b: fuse(list(val_truth[b]), val_seq_stream[b], val_str_stream[b],
                                        val_dom_stream[b], val_s_seq,
                                        tau, alpha, beta, gamma, tau_d, delta)
                                for b in BRANCHES
                            }
                            score = macro_fmax(cells, fused, val_truth)
                            if best is None or score > best[0]:
                                best = (score, tau, alpha, beta, gamma, tau_d, delta)
                                print("  new best: macro-Fmax={:.4f}  tau={} alpha={} beta={} "
                                      "gamma={} tau_d={} delta={}".format(*best))
    macro, tau, alpha, beta, gamma, tau_d, delta = best
    with open(r / "fusion_params.json", "w", encoding="utf-8") as f:
        json.dump({"tau": tau, "alpha": alpha, "beta": beta, "gamma": gamma,
                   "tau_d": tau_d, "delta": delta,
                   "val_macro_fmax": round(macro, 4)}, f, indent=2)
    print("Fitted: tau={} alpha={} beta={} gamma={} tau_d={} delta={} (val macro-Fmax {:.4f})".format(
        tau, alpha, beta, gamma, tau_d, delta, macro))

    # ---------------- test-side streams (full-train knowledge) -------------
    train_set = set(train)
    train_terms = {a: prop[a] for a in train if a in prop}
    test_seq_hits = load_clean_hits(d / "seq_hits.tsv", train_set)
    test_str_hits = load_clean_hits(d / "struct_hits.tsv", train_set)
    test_s_seq, _ = parse_raw_seq(d / "bin_hits.tsv", train_set, transfer_e)

    seq_stream = split_by_branch(transfer_stream(test_seq_hits, train_terms), dag)
    str_stream = split_by_branch(transfer_stream(test_str_hits, train_terms), dag)
    train_probs = domain_stats(train, prop, domains)
    dom_stream = split_by_branch(domain_stream(test, domains, train_probs, ip2go), dag)

    fused = {
        b: fuse(test, seq_stream[b], str_stream[b], dom_stream[b],
                test_s_seq, tau, alpha, beta, gamma, tau_d, delta)
        for b in BRANCHES
    }
    # ablation arm SD = sequence + domain, structure OFF (review point 2):
    # C minus SD isolates the incremental value of structure over seq+domain.
    empty_str = {b: {} for b in BRANCHES}
    seqdom = {
        b: fuse(test, seq_stream[b], empty_str[b], dom_stream[b],
                test_s_seq, tau, alpha, beta, gamma, tau_d, delta)
        for b in BRANCHES
    }
    nc = write_preds(r / "pred_armC.tsv", fused)
    ns = write_preds(r / "pred_armS.tsv", str_stream)
    nd = write_preds(r / "pred_armD.tsv", dom_stream)
    nsd = write_preds(r / "pred_armSD.tsv", seqdom)
    print("Arm C (fused): {} terms; S (struct): {}; D (domain): {}; "
          "SD (seq+domain): {}".format(nc, ns, nd, nsd))


if __name__ == "__main__":
    main()
