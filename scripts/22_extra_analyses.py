"""Step 22 (analysis) - naive baseline, per-species breakdown, family analysis.

A. Naive frequency baseline (CAFA convention): every term scored by its
   frequency among training proteins, identically for every test protein.
B. Per-species breakdown of test-set Fmax for the sequence arm and the
   gated fusion.
C. Family-level analysis: which InterPro entries account for the
   twilight-zone gains, the composition of the 246 sequence-orphan proteins,
   and dark-set versus annotated-set family composition.

Outputs:
  results/naive_baseline.csv
  results/species_breakdown.md
  results/family_analysis.md
"""
import csv
import sys
from collections import Counter, defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import data_dir, http_session, load_config, results_dir
from pis.eval import fmax, protein_curve
from pis.go import GoDag, parse_gaf
from pis.streams import BRANCHES, load_domains

BIN_ORDER = ["lt30", "30to50", "50to80", "ge80"]


def load_preds(path):
    preds = defaultdict(dict)
    with open(path, encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, term, score = line.rstrip("\n").split("\t")
            preds[acc][term] = float(score)
    return preds


def protein_best_f1(pred, truth):
    npred, ncorr = protein_curve(pred, truth)
    best = 0.0
    nt = len(truth)
    for j in range(len(npred)):
        if npred[j] == 0:
            continue
        pr = ncorr[j] / npred[j]
        rc = ncorr[j] / nt
        if pr + rc > 0:
            best = max(best, 2 * pr * rc / (pr + rc))
    return best


def ipr_name(session, ipr, cache):
    if ipr in cache:
        return cache[ipr]
    try:
        js = session.get("https://www.ebi.ac.uk/interpro/api/entry/interpro/"
                         + ipr, timeout=30).json()
        name = js["metadata"]["name"]
        if isinstance(name, dict):
            name = name.get("short") or name.get("name") or ipr
    except Exception:
        name = ipr
    cache[ipr] = name
    return name


def main():
    cfg = load_config()
    d = data_dir(cfg)
    r = results_dir(cfg)
    dag = GoDag(d / "go-basic.obo")
    session = http_session()

    train, test, bins = [], [], {}
    with open(d / "split.csv", encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, split, b = line.rstrip("\n").split(",")
            (train if split == "train" else test).append(acc)
            if split == "test":
                bins[acc] = b
    raw_ann = defaultdict(set)
    with open(d / "annotations.tsv", encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, term, _a = line.rstrip("\n").split("\t")
            raw_ann[acc].add(term)
    prop = {a: dag.propagate(ts) for a, ts in raw_ann.items()}

    truth = {b: {} for b in BRANCHES}
    for acc in test:
        for b in BRANCHES:
            sub = {t for t in prop.get(acc, ()) if dag.branch(t) == b}
            if sub:
                truth[b][acc] = sub

    # ---------------- A. naive frequency baseline --------------------------
    n_train = len(train)
    freq = {b: Counter() for b in BRANCHES}
    for a in train:
        for t in prop.get(a, ()):
            b = dag.branch(t)
            if b in freq:
                freq[b][t] += 1
    naive = {b: {t: c / n_train for t, c in freq[b].items()} for b in BRANCHES}

    rows = []
    for b in BRANCHES:
        for lbl in BIN_ORDER:
            prots = sorted(a for a in truth[b] if bins.get(a) == lbl)
            preds = {a: naive[b] for a in prots}
            f, t = fmax(prots, preds, truth[b])
            rows.append([b, lbl, len(prots), "{:.4f}".format(f)])
            print("naive {:22s} {:7s} n={:4d} Fmax={:.4f}".format(b, lbl, len(prots), f))
    with open(r / "naive_baseline.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["branch", "bin", "n", "fmax"])
        w.writerows(rows)

    # ---------------- B. per-species breakdown -----------------------------
    species_of = {}
    for sp in cfg["species"]:
        gaf = d / "gaf" / "goa_{}.gaf.gz".format(sp["name"])
        for acc, _t, _a in parse_gaf(gaf, cfg["evidence_codes"]):
            species_of.setdefault(acc, sp["name"])
    arm_a = load_preds(r / "pred_armA.tsv")
    arm_c = load_preds(r / "pred_armC.tsv")
    with open(r / "species_breakdown.md", "w", encoding="utf-8", newline="\n") as f:
        f.write("# Test-set Fmax by species (all identity bins pooled)\n\n")
        f.write("| species | branch | n | sequence (A) | gated fusion (C) |\n")
        f.write("|---|---|---|---|---|\n")
        for sp in cfg["species"]:
            for b in BRANCHES:
                prots = sorted(a for a in truth[b] if species_of.get(a) == sp["name"])
                pa = {a: {t: s for t, s in arm_a.get(a, {}).items()
                          if dag.branch(t) == b} for a in prots}
                pc = {a: {t: s for t, s in arm_c.get(a, {}).items()
                          if dag.branch(t) == b} for a in prots}
                fa, _ = fmax(prots, pa, truth[b])
                fc, _ = fmax(prots, pc, truth[b])
                f.write("| {} | {} | {} | {:.3f} | {:.3f} |\n".format(
                    sp["name"], b, len(prots), fa, fc))
                print("species {:6s} {:22s} n={:4d} A={:.3f} C={:.3f}".format(
                    sp["name"], b, len(prots), fa, fc))

    # ---------------- C. family analysis -----------------------------------
    domains = load_domains(d / "domains.tsv")
    dark_domains = load_domains(d / "dark" / "domains.tsv")
    name_cache = {}

    # C1: per-protein fusion gain in the <30% bins, grouped by InterPro entry
    gain_by_ipr = defaultdict(list)
    for b in BRANCHES:
        for acc in truth[b]:
            if bins.get(acc) != "lt30":
                continue
            pa = {t: s for t, s in arm_a.get(acc, {}).items() if dag.branch(t) == b}
            pc = {t: s for t, s in arm_c.get(acc, {}).items() if dag.branch(t) == b}
            delta = protein_best_f1(pc, truth[b][acc]) - protein_best_f1(pa, truth[b][acc])
            for ipr in domains.get(acc, ()):
                gain_by_ipr[ipr].append(delta)
    ranked_gain = sorted(((ipr, len(v), sum(v) / len(v))
                          for ipr, v in gain_by_ipr.items() if len(v) >= 5),
                         key=lambda x: -x[2])[:10]

    # C2: composition of the 246 sequence-orphan (rescued) proteins
    rescued = sorted(set(load_preds(r / "pred_armB.tsv")) - set(arm_a))
    resc_counter = Counter(ipr for a in rescued for ipr in domains.get(a, ()))
    resc_nodomain = sum(1 for a in rescued if not domains.get(a))

    # C3: dark vs annotated family composition
    kb = set((d / "final_accessions.txt").read_text(encoding="utf-8").split())
    dark = set((d / "dark" / "accessions.txt").read_text(encoding="utf-8").split())
    kb_counter = Counter(ipr for a in kb for ipr in domains.get(a, ()))
    dark_counter = Counter(ipr for a in dark for ipr in dark_domains.get(a, ()))
    dark_nodomain = sum(1 for a in dark if not dark_domains.get(a))

    with open(r / "family_analysis.md", "w", encoding="utf-8", newline="\n") as f:
        f.write("# Family-level analysis\n\n")
        f.write("## InterPro entries with the largest mean fusion gain, <30% identity bins (n >= 5 proteins)\n\n")
        f.write("| InterPro | name | proteins | mean per-protein F1 gain (C - A) |\n|---|---|---|---|\n")
        for ipr, n, g in ranked_gain:
            f.write("| {} | {} | {} | {:+.3f} |\n".format(ipr, ipr_name(session, ipr, name_cache), n, g))
        f.write("\n## Composition of the {} proteins with no sequence neighbours\n\n".format(len(rescued)))
        f.write("{} of {} have no InterPro entry at all.\n\n".format(resc_nodomain, len(rescued)))
        f.write("| InterPro | name | proteins |\n|---|---|---|\n")
        for ipr, n in resc_counter.most_common(10):
            f.write("| {} | {} | {} |\n".format(ipr, ipr_name(session, ipr, name_cache), n))
        f.write("\n## Most common InterPro entries: annotated set vs dark set\n\n")
        f.write("{} of {} dark proteins have no InterPro entry.\n\n".format(dark_nodomain, len(dark)))
        f.write("| rank | annotated set | n | dark set | n |\n|---|---|---|---|---|\n")
        top_kb = kb_counter.most_common(10)
        top_dark = dark_counter.most_common(10)
        for i in range(10):
            ka, kn = top_kb[i] if i < len(top_kb) else ("-", "")
            da, dn = top_dark[i] if i < len(top_dark) else ("-", "")
            f.write("| {} | {} ({}) | {} | {} ({}) | {} |\n".format(
                i + 1, ka, ipr_name(session, ka, name_cache) if ka != "-" else "-",
                kn, da, ipr_name(session, da, name_cache) if da != "-" else "-", dn))
    print("rescued no-domain: {}/{}; dark no-domain: {}/{}".format(
        resc_nodomain, len(rescued), dark_nodomain, len(dark)))
    print("Wrote results/naive_baseline.csv, species_breakdown.md, family_analysis.md")


if __name__ == "__main__":
    main()
