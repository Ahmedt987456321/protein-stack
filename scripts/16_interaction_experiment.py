"""Step 16 (v0.4) - does an interaction stream add signal on top of fusion?

Arm I: guilt-by-association - transfer GO terms from STRING partners
(combined score / 1000 as similarity). Arm C+I: noisy-OR of the v0.2 fused
predictions with the interaction stream, weight w_i fitted on the SAME
validation split used in v0.2 (partners restricted to subtrain; never test).

Gate (v0.4): macro-Fmax(C+I) >= macro-Fmax(C), and no statistically
significant dilution in any branch x bin cell (per-cell bootstrap, as v0.2).

Outputs:
  results/pred_armCI.tsv, pred_armI.tsv
  results/report_v04.md   (also folds in STRING + pocket summaries)
"""
import json
import random
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import data_dir, load_config, results_dir
from pis.eval import fmax, fmax_from_curves, protein_curve
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

GRID_WI = [0.2, 0.4, 0.6, 0.8, 1.0]
MIN_CELL = 10
N_BOOT = 1000


def load_string(path):
    adj = defaultdict(list)
    with open(path, encoding="utf-8") as f:
        next(f)
        for line in f:
            a, b, score = line.rstrip("\n").split("\t")
            s = int(score) / 1000.0
            adj[a].append((b, s))
            adj[b].append((a, s))
    return adj


def interaction_hits(adj, queries, allowed_targets):
    hits = {}
    for q in queries:
        partners = [(t, s) for t, s in adj.get(q, ()) if t in allowed_targets]
        if partners:
            hits[q] = partners
    return hits


def or_combine(base, extra, w):
    """Per-protein noisy-OR of two term->score dicts."""
    out = dict(base)
    for t, s in extra.items():
        prev = out.get(t, 0.0)
        out[t] = 1.0 - (1.0 - prev) * (1.0 - w * s)
    return out


def main():
    cfg = load_config()
    d = data_dir(cfg)
    r = results_dir(cfg)
    edges_cfg = cfg["split"]["bin_edges"]
    transfer_e = cfg["search"]["evalue"]

    dag = GoDag(d / "go-basic.obo")
    params = json.load(open(r / "fusion_params.json", encoding="utf-8"))
    adj = load_string(d / "string_edges.tsv")
    print("STRING adjacency: {} proteins with >=1 partner".format(len(adj)))

    train, test = [], []
    bins = {}
    with open(d / "split.csv", encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, split, b = line.rstrip("\n").split(",")
            if split == "train":
                train.append(acc)
            else:
                test.append(acc)
                bins[acc] = b
    train, test = sorted(train), sorted(test)

    raw_ann = defaultdict(set)
    with open(d / "annotations.tsv", encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, term, _a = line.rstrip("\n").split("\t")
            raw_ann[acc].add(term)
    prop = {acc: dag.propagate(terms) for acc, terms in raw_ann.items()}

    domains = load_domains(d / "domains.tsv")
    ip2go = load_interpro2go(d / "interpro2go.txt", dag)

    # ---------------- fit w_i on the v0.2 validation split ------------------
    rng = random.Random(cfg["seed"] + 1)          # SAME split as script 09
    shuffled = list(train)
    rng.shuffle(shuffled)
    n_val = round(0.2 * len(shuffled))
    val = sorted(shuffled[:n_val])
    subtrain = sorted(shuffled[n_val:])
    sub_set = set(subtrain)
    sub_terms = {a: prop[a] for a in subtrain if a in prop}

    val_s_seq, val_seq_transfer = parse_raw_seq(d / "val_seq_hits_raw.tsv", sub_set, transfer_e)
    val_struct = parse_raw_struct(d / "val_struct_hits_raw.tsv", sub_set)
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
    val_bins = {acc: bin_of(val_s_seq.get(acc, 0.0), edges_cfg) for acc in val}
    val_cells = []
    for b in BRANCHES:
        for lbl in BIN_LABELS:
            members = sorted(a for a in val_truth[b] if val_bins[a] == lbl)
            if len(members) >= MIN_CELL:
                val_cells.append((b, members))

    val_c = {
        b: fuse(list(val_truth[b]), val_seq_stream[b], val_str_stream[b],
                val_dom_stream[b], val_s_seq,
                params["tau"], params["alpha"], params["beta"],
                params["gamma"], params["tau_d"], params["delta"])
        for b in BRANCHES
    }
    val_i_hits = interaction_hits(adj, val, sub_set)
    val_i_stream = split_by_branch(transfer_stream(val_i_hits, sub_terms), dag)
    print("Validation proteins with subtrain STRING partners: {}".format(len(val_i_hits)))

    def val_macro(preds_by_branch):
        scores = []
        for b, members in val_cells:
            f, _ = fmax(members, preds_by_branch[b], val_truth[b])
            scores.append(f)
        return sum(scores) / len(scores)

    base_macro = val_macro(val_c)
    best = (base_macro, 0.0)
    for wi in GRID_WI:
        fused = {
            b: {acc: or_combine(val_c[b].get(acc, {}), val_i_stream[b].get(acc, {}), wi)
                for acc in val_truth[b]}
            for b in BRANCHES
        }
        m = val_macro(fused)
        print("  w_i={}: val macro-Fmax {:.4f} (C alone {:.4f})".format(wi, m, base_macro))
        if m > best[0]:
            best = (m, wi)
    w_i = best[1]
    print("Fitted w_i = {} (val macro {:.4f})".format(w_i, best[0]))

    # ---------------- test-side arms ----------------------------------------
    test_truth = {b: {} for b in BRANCHES}
    for acc in test:
        for b in BRANCHES:
            sub = {t for t in prop.get(acc, ()) if dag.branch(t) == b}
            if sub:
                test_truth[b][acc] = sub

    train_set = set(train)
    train_terms = {a: prop[a] for a in train if a in prop}
    test_i_hits = interaction_hits(adj, test, train_set)
    test_i_stream = split_by_branch(transfer_stream(test_i_hits, train_terms), dag)
    print("Test proteins with train STRING partners: {}".format(len(test_i_hits)))

    pred_c = defaultdict(dict)
    with open(r / "pred_armC.tsv", encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, term, score = line.rstrip("\n").split("\t")
            pred_c[acc][term] = float(score)
    pred_c_branch = {
        b: {acc: {t: s for t, s in pred_c.get(acc, {}).items() if dag.branch(t) == b}
            for acc in test_truth[b]}
        for b in BRANCHES
    }
    pred_ci_branch = {
        b: {acc: or_combine(pred_c_branch[b].get(acc, {}),
                            test_i_stream[b].get(acc, {}), w_i)
            for acc in test_truth[b]}
        for b in BRANCHES
    }

    def write_arm(path, by_branch):
        merged = defaultdict(dict)
        for per in by_branch.values():
            for acc, scores in per.items():
                merged[acc].update(scores)
        with open(path, "w", encoding="utf-8", newline="\n") as f:
            f.write("accession\tgo_term\tscore\n")
            for acc in sorted(merged):
                for term, score in sorted(merged[acc].items()):
                    f.write("{}\t{}\t{:.4f}\n".format(acc, term, score))
    write_arm(r / "pred_armCI.tsv", pred_ci_branch)
    write_arm(r / "pred_armI.tsv", {b: test_i_stream[b] for b in BRANCHES})

    # ---------------- evaluate + gate ---------------------------------------
    table = {}
    arms = {"C": pred_c_branch, "I": {b: test_i_stream[b] for b in BRANCHES}, "CI": pred_ci_branch}
    cells = []
    for b in BRANCHES:
        for lbl in BIN_LABELS:
            members = sorted(a for a in test_truth[b] if bins.get(a) == lbl)
            if members:
                cells.append((b, lbl, members))
            for arm, preds in arms.items():
                f, _ = fmax(members, preds[b], test_truth[b])
                table[(arm, b, lbl)] = (f, len(members))

    macro = {arm: sum(table[(arm, b, lbl)][0] for b, lbl, _m in cells) / len(cells)
             for arm in arms}

    rng2 = random.Random(cfg["seed"] + 5)

    def cell_ci(b, members, arm_hi, arm_lo):
        curves = {
            arm: [protein_curve(arms[arm][b].get(acc, {}), test_truth[b][acc]) for acc in members]
            for arm in (arm_hi, arm_lo)
        }
        sizes = [len(test_truth[b][acc]) for acc in members]
        deltas = []
        for _ in range(N_BOOT):
            idx = [rng2.randrange(len(members)) for _ in members]
            fh, _ = fmax_from_curves([curves[arm_hi][i] for i in idx], [sizes[i] for i in idx])
            fl, _ = fmax_from_curves([curves[arm_lo][i] for i in idx], [sizes[i] for i in idx])
            deltas.append(fh - fl)
        deltas.sort()
        return deltas[int(0.025 * N_BOOT)], deltas[int(0.975 * N_BOOT) - 1]

    dilution_checks = []
    improvements = []
    for b, lbl, members in cells:
        delta_pt = table[("CI", b, lbl)][0] - table[("C", b, lbl)][0]
        if delta_pt < 0:
            lo, hi = cell_ci(b, members, "CI", "C")
            dilution_checks.append((b, lbl, delta_pt, lo, hi, hi < 0))
        elif delta_pt > 0.005:
            improvements.append((b, lbl, delta_pt))
    diluted = [(b, lbl) for b, lbl, _d, _l, _h, sig in dilution_checks if sig]
    g_macro = macro["CI"] >= macro["C"]
    g_no_dilution = not diluted
    verdict = "PASSED" if (g_macro and g_no_dilution) else "NOT PASSED"

    # ---------------- report -------------------------------------------------
    n_string = sum(1 for _ in open(d / "string_edges.tsv", encoding="utf-8")) - 1
    pockets_line = ""
    pk = d / "dark" / "pockets.tsv"
    if pk.exists():
        rows = [line.split("\t") for line in open(pk, encoding="utf-8")]
        n_drug = sum(1 for c in rows if len(c) >= 3 and float(c[2]) >= 0.5)
        pockets_line = ("fpocket ran over {} dark-protein models; {} have a pocket with "
                        "druggability >= 0.5 (see results/druggability_top.md).").format(
                            len(rows), n_drug)

    with open(r / "report_v04.md", "w", encoding="utf-8", newline="\n") as f:
        f.write("# v0.4 - interactions + pockets - report\n\n")
        f.write("STRING edges (combined score >= {}): {} among our {} proteins. "
                "Interaction arm I = guilt-by-association transfer from STRING "
                "partners; C+I = noisy-OR on top of the v0.2 fused arm with "
                "w_i = {} fitted on the v0.2 validation split.\n\n".format(
                    cfg["string"]["score_min"], n_string,
                    len(set((d / "final_accessions.txt").read_text().split())), w_i))
        if pockets_line:
            f.write(pockets_line + "\n\n")
        f.write("## Gate\n\n")
        f.write("- macro-Fmax: C {:.4f} -> C+I {:.4f} ({}); I alone {:.4f}\n".format(
            macro["C"], macro["CI"], "PASS" if g_macro else "FAIL", macro["I"]))
        f.write("- dilution: {}\n".format(
            "PASS" if g_no_dilution else "FAIL: " + ", ".join("{}/{}".format(b, l) for b, l in diluted)))
        for b, lbl, dpt, lo, hi, sig in dilution_checks:
            f.write("  - {}/{}: point {:+.4f}, CI [{:+.4f}, {:+.4f}] -> {}\n".format(
                b, lbl, dpt, lo, hi, "SIGNIFICANT" if sig else "not significant"))
        f.write("\n**Phase 0.4 gate: {}**\n\n".format(verdict))
        f.write("## Fmax by branch / bin\n\n")
        f.write("| branch | bin | n | I | C | C+I | delta |\n|---|---|---:|---:|---:|---:|---:|\n")
        for b in BRANCHES:
            for lbl in BIN_LABELS:
                fi = table[("I", b, lbl)][0]
                fc = table[("C", b, lbl)][0]
                fci = table[("CI", b, lbl)][0]
                f.write("| {} | {} | {} | {:.4f} | {:.4f} | {:.4f} | {:+.4f} |\n".format(
                    b, lbl, table[("C", b, lbl)][1], fi, fc, fci, fci - fc))
        if improvements:
            f.write("\nCells improved by > 0.005: " + ", ".join(
                "{}/{} ({:+.4f})".format(b, l, dv) for b, l, dv in improvements) + "\n")
        f.write("\nEvidence type of interaction transfer: COMPUTATIONAL (STRING "
                "combined scores include predicted channels).\n")

    for b in BRANCHES:
        for lbl in BIN_LABELS:
            print("{:22s} {:7s} n={:5d}  I={:.4f}  C={:.4f}  C+I={:.4f}".format(
                b, lbl, table[("C", b, lbl)][1], table[("I", b, lbl)][0],
                table[("C", b, lbl)][0], table[("CI", b, lbl)][0]))
    print()
    print("macro: I={:.4f}  C={:.4f}  C+I={:.4f}".format(macro["I"], macro["C"], macro["CI"]))
    print("PHASE 0.4 GATE: " + verdict)
    print("Wrote results/report_v04.md")


if __name__ == "__main__":
    main()
