"""Step 12 (v0.3) - dark-proteome sweep: generate ranked function hypotheses
for proteins with no experimentally supported annotation.

Knowledge base = all 9K annotated proteins (train+test; no held-out hygiene is
needed because dark proteins have no experimental truth to leak). Streams and
the fitted fusion parameters come from v0.2 unchanged.

Corroboration: dark proteins' electronic (IEA-class) annotations were held
aside in step 11 and are NEVER inputs to prediction. Agreement between our
top hypothesis and a protein's electronic annotations is corroboration, not
proof - both are computational - but a HIGH-confidence tier that agrees far
more often than a shuffled baseline is behaving like a signal, not noise.

Gate (v0.3): HIGH-tier top-1 corroboration >= 2x the shuffled baseline, and
at least 100 HIGH-tier hypotheses produced.

Outputs:
  results/dark_hypotheses.tsv
  results/report_v03.md
"""
import json
import random
import shutil
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import data_dir, load_config, results_dir, run_tool
from pis.go import GoDag
from pis.streams import (
    BRANCHES,
    domain_stats,
    domain_stream,
    fuse,
    load_domains,
    load_interpro2go,
    parse_raw_seq,
    parse_raw_struct,
    specific_terms,
    split_by_branch,
    transfer_stream,
)

TIER_HIGH = 0.7
TIER_MED = 0.45
N_SHUFFLE = 100


def tier_of(score):
    if score >= TIER_HIGH:
        return "HIGH"
    if score >= TIER_MED:
        return "MEDIUM"
    return "LOW"


def ensure_dark_searches(cfg, d, dd):
    th = str(cfg["tools"]["threads"])
    if not (dd / "seq_hits_raw.tsv").exists():
        run_tool(cfg, "mmseqs", [
            "easy-search", dd / "dark.fasta", d / "final.fasta",
            dd / "seq_hits_raw.tsv", dd / "tmp_seq",
            "--format-output", "query,target,pident,evalue",
            "-e", "10", "-s", "7.5",
            "--max-seqs", str(cfg["search"]["max_seqs"]), "--threads", th,
        ])
        shutil.rmtree(dd / "tmp_seq", ignore_errors=True)
    if not (dd / "struct_hits_raw.tsv").exists():
        run_tool(cfg, "foldseek", [
            "easy-search", dd / "structures", d / "structures",
            dd / "struct_hits_raw.tsv", dd / "tmp_struct",
            "--format-output", "query,target,alntmscore,bits,evalue",
            "-e", str(cfg["search"]["evalue"]),
            "--max-seqs", str(cfg["search"]["max_seqs"]), "--threads", th,
        ])
        shutil.rmtree(dd / "tmp_struct", ignore_errors=True)


def main():
    cfg = load_config()
    d = data_dir(cfg)
    dd = d / "dark"
    r = results_dir(cfg)
    transfer_e = cfg["search"]["evalue"]
    min_score = cfg["dark"]["min_score"]
    top_n = cfg["dark"]["top_per_branch"]

    dag = GoDag(d / "go-basic.obo")
    params = json.load(open(r / "fusion_params.json", encoding="utf-8"))
    print("Fusion params: {}".format(params))

    kb = (d / "final_accessions.txt").read_text(encoding="utf-8").split()
    kb_set = set(kb)
    dark = (dd / "accessions.txt").read_text(encoding="utf-8").split()
    print("Dark proteins: {}; knowledge base: {}".format(len(dark), len(kb)))

    ensure_dark_searches(cfg, d, dd)

    # ---- knowledge-base annotations (propagated) --------------------------
    raw_ann = defaultdict(set)
    with open(d / "annotations.tsv", encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, term, _aspect = line.rstrip("\n").split("\t")
            if acc in kb_set:
                raw_ann[acc].add(term)
    prop = {acc: dag.propagate(terms) for acc, terms in raw_ann.items()}

    # ---- streams -----------------------------------------------------------
    s_seq, seq_transfer = parse_raw_seq(dd / "seq_hits_raw.tsv", kb_set, transfer_e)
    struct_hits = parse_raw_struct(dd / "struct_hits_raw.tsv", kb_set)
    kb_domains = load_domains(d / "domains.tsv")
    dark_domains = load_domains(dd / "domains.tsv")
    ip2go = load_interpro2go(d / "interpro2go.txt", dag)

    seq_stream = split_by_branch(transfer_stream(seq_transfer, prop), dag)
    str_stream = split_by_branch(transfer_stream(struct_hits, prop), dag)
    kb_probs = domain_stats(kb, prop, kb_domains)
    dom_stream = split_by_branch(domain_stream(dark, dark_domains, kb_probs, ip2go), dag)

    fused = {
        b: fuse(dark, seq_stream[b], str_stream[b], dom_stream[b], s_seq,
                params["tau"], params["alpha"], params["beta"],
                params["gamma"], params["tau_d"], params["delta"])
        for b in BRANCHES
    }

    # provenance: which streams touch each protein at all
    has_seq = set(seq_transfer)
    has_struct = set(struct_hits)
    has_dom = {a for a in dark if dark_domains.get(a)}
    structure_only = [a for a in dark
                      if a in has_struct and a not in has_seq and a not in has_dom]

    # ---- hypothesis feed ---------------------------------------------------
    n_rows = 0
    tier_counts = defaultdict(int)
    top_mf = {}  # acc -> (term, score) best specific MF hypothesis
    with open(r / "dark_hypotheses.tsv", "w", encoding="utf-8", newline="\n") as f:
        f.write("accession\tbranch\tgo_term\tscore\ttier\tstreams\tstructure_only\n")
        for b in BRANCHES:
            for acc in sorted(fused[b]):
                picks = specific_terms(fused[b][acc], dag, min_score, top_n)
                for term, score in picks:
                    streams = []
                    if seq_stream[b].get(acc, {}).get(term):
                        streams.append("seq")
                    if str_stream[b].get(acc, {}).get(term):
                        streams.append("struct")
                    if dom_stream[b].get(acc, {}).get(term):
                        streams.append("dom")
                    t = tier_of(score)
                    tier_counts[t] += 1
                    f.write("{}\t{}\t{}\t{:.4f}\t{}\t{}\t{}\n".format(
                        acc, b, term, score, t, "+".join(streams),
                        int(acc in structure_only)))
                    n_rows += 1
                if b == "molecular_function" and picks:
                    top_mf[acc] = picks[0]
    print("Hypotheses written: {} ({})".format(
        n_rows, ", ".join("{}={}".format(k, tier_counts[k]) for k in ("HIGH", "MEDIUM", "LOW"))))
    print("Structure-only proteins (no seq hits, no domains): {}".format(len(structure_only)))

    # ---- corroboration against held-aside electronic annotations ----------
    iea = defaultdict(set)
    with open(dd / "iea_annotations.tsv", encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, term, _aspect = line.rstrip("\n").split("\t")
            iea[acc].add(term)
    iea_prop = {}
    for acc in dark:
        terms = dag.propagate(iea.get(acc, set()))
        mf = {t for t in terms if dag.branch(t) == "molecular_function"}
        if mf:
            iea_prop[acc] = mf

    checkable = {acc: tv for acc, tv in top_mf.items() if acc in iea_prop}
    agree = {t: [0, 0] for t in ("HIGH", "MEDIUM", "LOW")}  # [agree, total]
    for acc, (term, score) in checkable.items():
        t = tier_of(score)
        agree[t][1] += 1
        if term in iea_prop[acc]:
            agree[t][0] += 1

    rng = random.Random(cfg["seed"] + 4)
    accs_c = sorted(checkable)
    terms_c = [checkable[a][0] for a in accs_c]
    base_hits = 0
    for _ in range(N_SHUFFLE):
        perm = terms_c[:]
        rng.shuffle(perm)
        base_hits += sum(1 for a, t in zip(accs_c, perm) if t in iea_prop[a])
    baseline = base_hits / (N_SHUFFLE * len(accs_c)) if accs_c else 0.0

    high_rate = agree["HIGH"][0] / agree["HIGH"][1] if agree["HIGH"][1] else 0.0
    g_ratio = (high_rate / baseline) if baseline > 0 else float("inf")
    g_pass = g_ratio >= 2.0 and tier_counts["HIGH"] >= 100
    verdict = "PASSED" if g_pass else "NOT PASSED"

    # ---- report ------------------------------------------------------------
    with open(r / "report_v03.md", "w", encoding="utf-8", newline="\n") as f:
        f.write("# v0.3 - dark-proteome sweep - report\n\n")
        f.write("Dark set: **{} proteins** (GOA proteins with zero experimental GO "
                "evidence), knowledge base: {} annotated proteins, fusion params "
                "from v0.2: `{}`.\n\n".format(len(dark), len(kb), params))
        f.write("## Hypothesis feed\n\n")
        f.write("- {} hypotheses across {} proteins "
                "(HIGH {}, MEDIUM {}, LOW {})\n".format(
                    n_rows, len({a for b in BRANCHES for a in fused[b]}),
                    tier_counts["HIGH"], tier_counts["MEDIUM"], tier_counts["LOW"]))
        f.write("- **Structure-only proteins: {}** - no sequence neighbours, no "
                "InterPro domains; AlphaFold+Foldseek is the only evidence. "
                "These are the highest-novelty targets.\n\n".format(len(structure_only)))
        f.write("## Corroboration vs held-aside electronic annotations "
                "(top MF hypothesis, n={} checkable)\n\n".format(len(checkable)))
        f.write("| tier | agree | total | rate |\n|---|---:|---:|---:|\n")
        for t in ("HIGH", "MEDIUM", "LOW"):
            a, n = agree[t]
            f.write("| {} | {} | {} | {} |\n".format(
                t, a, n, "{:.1%}".format(a / n) if n else "-"))
        f.write("| shuffled baseline | - | - | {:.1%} |\n\n".format(baseline))
        f.write("Electronic annotations are themselves computational; agreement is "
                "corroboration, not validation. Disagreements in the HIGH tier are "
                "the interesting review queue.\n\n")
        f.write("**Gate** (HIGH-tier corroboration >= 2x shuffled baseline, and "
                ">= 100 HIGH hypotheses): HIGH {:.1%} vs baseline {:.1%} "
                "(ratio {:.1f}x), HIGH count {} -> **{}**\n".format(
                    high_rate, baseline, g_ratio, tier_counts["HIGH"], verdict))
    print("Corroboration: " + ", ".join(
        "{} {}/{}".format(t, agree[t][0], agree[t][1]) for t in ("HIGH", "MEDIUM", "LOW")))
    print("HIGH rate {:.1%} vs shuffled baseline {:.1%} (ratio {:.1f}x)".format(
        high_rate, baseline, g_ratio))
    print("PHASE 0.3 GATE: " + verdict)
    print("Wrote results/dark_hypotheses.tsv and results/report_v03.md")


if __name__ == "__main__":
    main()
