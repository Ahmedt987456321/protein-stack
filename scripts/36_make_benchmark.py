"""Exploration B1 - generate the timesplit-go benchmark package.

Writes benchmark/cohorts/, truth/, ic/, baselines/ (SS example
submissions), baselines.md, and metadata.json from the cached per-horizon
corpora. Truth is the propagated current experimental annotation set per
cohort protein, frozen at generation time. IC tables (information content,
-log2 of term frequency in the horizon's past corpus) are frozen per
horizon and drive the primary metric, mean information gain.
"""
import csv
import datetime
import json
import math
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import data_dir, load_config
from pis.go import GoDag, parse_gaf
from pis.streams import (
    BRANCHES,
    fuse,
    parse_raw_seq,
    parse_raw_struct,
    specific_terms,
    split_by_branch,
    transfer_stream,
)

HORIZONS = {"2019": "timesplit_h2019", "2021": "timesplit_h2021",
            "2023": "timesplit_h2023"}


def main():
    cfg = load_config()
    d = data_dir(cfg)
    bdir = Path("benchmark")
    for sub in ("cohorts", "truth", "ic", "baselines"):
        (bdir / sub).mkdir(parents=True, exist_ok=True)
    dag = GoDag(d / "go-basic.obo")
    params = json.load(open(Path("results") / "fusion_params.json",
                            encoding="utf-8"))

    kb = set((d / "final_accessions.txt").read_text(encoding="utf-8").split())
    dark = set((d / "dark" / "accessions.txt").read_text(encoding="utf-8").split())
    universe = kb | dark

    current_exp = defaultdict(set)
    for sp in cfg["species"]:
        gaf = d / "gaf" / "goa_{}.gaf.gz".format(sp["name"])
        for acc, term, _a in parse_gaf(gaf, cfg["evidence_codes"]):
            current_exp[acc].add(term)

    meta = {"generated": datetime.date.today().isoformat(),
            "go_obo": cfg["go_obo_url"], "horizons": {}}
    baseline_rows = []
    for year, tsdir in HORIZONS.items():
        ts = d / tsdir
        past_exp = defaultdict(set)
        releases = []
        for gaf in sorted(ts.glob("goa_*.gaf.*.gz")):
            releases.append(gaf.name)
            for acc, term, _a in parse_gaf(gaf, cfg["evidence_codes"]):
                past_exp[acc].add(term)
        cohort = sorted(a for a in universe
                        if current_exp.get(a) and not past_exp.get(a))
        kb_past = sorted(a for a in kb - set(cohort) if past_exp.get(a))
        meta["horizons"][year] = {"past_releases": releases,
                                  "cohort_size": len(cohort)}

        (bdir / "cohorts" / "cohort_{}.txt".format(year)).write_text(
            "\n".join(cohort) + "\n", encoding="utf-8")
        truth = {b: {} for b in BRANCHES}
        with open(bdir / "truth" / "truth_{}.tsv".format(year), "w",
                  encoding="utf-8", newline="\n") as f:
            f.write("accession\tbranch\tgo_term\n")
            for acc in cohort:
                prop = dag.propagate(current_exp[acc])
                for b in BRANCHES:
                    sub = sorted(x for x in prop if dag.branch(x) == b)
                    if sub:
                        truth[b][acc] = set(sub)
                    for t in sub:
                        f.write("{}\t{}\t{}\n".format(acc, b, t))

        # IC from the past corpus (annotated non-cohort proteins)
        freq = Counter()
        n_prot = 0
        for acc in kb_past:
            terms = dag.propagate(past_exp[acc])
            if terms:
                n_prot += 1
                freq.update(terms)
        ic = {t: -math.log2(c / n_prot) for t, c in freq.items()}
        ic_max = max(ic.values())
        with open(bdir / "ic" / "ic_{}.tsv".format(year), "w",
                  encoding="utf-8", newline="\n") as f:
            f.write("go_term\tic_bits\n")
            for t, v in sorted(ic.items()):
                f.write("{}\t{:.4f}\n".format(t, v))

        # SS baseline predictions (same protocol as the manuscript)
        kb_past_set = set(kb_past)
        past_terms = {a: dag.propagate(past_exp[a]) for a in kb_past}
        s_seq, seq_tr = parse_raw_seq(ts / "seq_hits_raw.tsv", kb_past_set,
                                      cfg["search"]["evalue"])
        struct_hits = parse_raw_struct(ts / "struct_hits_raw.tsv", kb_past_set)
        seq_stream = split_by_branch(transfer_stream(seq_tr, past_terms), dag)
        str_stream = split_by_branch(transfer_stream(struct_hits, past_terms), dag)
        empty = {b: {} for b in BRANCHES}
        fused = {b: fuse(cohort, seq_stream[b], str_stream[b], empty[b], s_seq,
                         params["tau"], params["alpha"], params["beta"],
                         params["gamma"], params["tau_d"], params["delta"])
                 for b in BRANCHES}
        ss_pred = {b: {} for b in BRANCHES}
        with open(bdir / "baselines" / "ss_{}.tsv".format(year), "w",
                  encoding="utf-8", newline="\n") as f:
            f.write("accession\tbranch\tgo_term\n")
            for b in BRANCHES:
                for acc in cohort:
                    picks = specific_terms(fused[b].get(acc, {}), dag,
                                           cfg["timesplit"]["min_confidence"], 1)
                    if picks:
                        ss_pred[b][acc] = picks[0][0]
                        f.write("{}\t{}\t{}\n".format(acc, b, picks[0][0]))

        # score both baselines with both metrics
        prior_term = {b: max(((t, c) for t, c in freq.items()
                              if dag.branch(t) == b), key=lambda x: x[1])[0]
                      for b in BRANCHES}
        for b in BRANCHES:
            scored = sorted(a for a in truth[b] if a in ss_pred[b])
            hits = [a for a in scored if ss_pred[b][a] in truth[b][a]]
            n_sc = len(scored)
            ss_p = len(hits) / n_sc if n_sc else 0.0
            ss_gain = (sum(ic.get(ss_pred[b][a], ic_max) for a in hits) / n_sc
                       if n_sc else 0.0)
            pt = prior_term[b]
            cts = [a for a in truth[b]]
            ph = [a for a in cts if pt in truth[b][a]]
            pr_p = len(ph) / len(cts)
            pr_gain = len(ph) * ic.get(pt, ic_max) / len(cts)
            baseline_rows.append([year, b, len(cohort), len(scored),
                                  "{:.3f}".format(ss_p), "{:.3f}".format(ss_gain),
                                  "{} ({})".format(pt, dag.name(pt)),
                                  "{:.3f}".format(pr_p), "{:.3f}".format(pr_gain)])
        print(year, "cohort", len(cohort), releases)

    with open(bdir / "baselines.md", "w", encoding="utf-8", newline="\n") as f:
        f.write("# Baseline results\n\n")
        f.write("Primary metric: mean information gain (bits) = IC of the "
                "predicted term when correct, 0 when wrong, averaged over "
                "scored proteins; IC tables frozen per horizon from the past "
                "corpus. Raw top-1 precision shown for interpretability - "
                "note the prior's high raw precision but near-zero "
                "information gain, which is why raw precision is not the "
                "ranking metric.\n\n")
        f.write("| horizon | branch | cohort | SS n | SS top-1 | SS gain | "
                "prior term | prior top-1 | prior gain |\n")
        f.write("|---|---|---|---|---|---|---|---|---|\n")
        for r in baseline_rows:
            f.write("| " + " | ".join(str(x) for x in r) + " |\n")

    (bdir / "metadata.json").write_text(json.dumps(meta, indent=2),
                                        encoding="utf-8")
    print("Wrote benchmark/ package")


if __name__ == "__main__":
    main()
