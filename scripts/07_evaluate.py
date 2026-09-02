"""Step 07 - CAFA-style Fmax per identity bin and GO branch, plus a bootstrap
confidence interval on the primary endpoint:

  delta Fmax (Arm B - Arm A) in the lt30 bin, molecular_function.

Outputs:
  results/summary.csv   one row per (arm, branch, bin)
  results/report.md     human-readable report with the verdict
"""
import csv
import random
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import data_dir, load_config, results_dir
from pis.eval import fmax, fmax_from_curves, protein_curve
from pis.go import GoDag

BRANCHES = ["molecular_function", "biological_process", "cellular_component"]
BIN_ORDER = ["lt30", "30to50", "50to80", "ge80"]
N_BOOT = 1000


def load_test_bins(path: Path):
    bins = {}
    with open(path, encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, split, b = line.rstrip("\n").split(",")
            if split == "test":
                bins[acc] = b
    return bins


def load_preds(path: Path):
    preds = defaultdict(dict)
    with open(path, encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, term, score = line.rstrip("\n").split("\t")
            preds[acc][term] = float(score)
    return preds


def main():
    cfg = load_config()
    d = data_dir(cfg)
    r = results_dir(cfg)

    dag = GoDag(d / "go-basic.obo")
    bins = load_test_bins(d / "split.csv")

    raw_ann = defaultdict(set)
    with open(d / "annotations.tsv", encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, term, _aspect = line.rstrip("\n").split("\t")
            if acc in bins:
                raw_ann[acc].add(term)

    # propagated truth per branch
    truth = {b: {} for b in BRANCHES}
    for acc, terms in raw_ann.items():
        prop = dag.propagate(terms)
        for b in BRANCHES:
            sub = {t for t in prop if dag.branch(t) == b}
            if sub:
                truth[b][acc] = sub

    arms = {"A": load_preds(r / "pred_armA.tsv"), "B": load_preds(r / "pred_armB.tsv")}

    rows = []
    fmax_table = {}  # (arm, branch, bin) -> fmax
    for branch in BRANCHES:
        for b in BIN_ORDER:
            proteins = sorted(acc for acc in truth[branch] if bins.get(acc) == b)
            for arm, preds in arms.items():
                branch_preds = {
                    acc: {t: s for t, s in preds.get(acc, {}).items() if dag.branch(t) == branch}
                    for acc in proteins
                }
                f, t = fmax(proteins, branch_preds, truth[branch])
                fmax_table[(arm, branch, b)] = f
                rows.append([arm, branch, b, len(proteins), "{:.4f}".format(f), "{:.2f}".format(t)])
                print("arm {}  {:22s} {:7s} n={:5d}  Fmax={:.4f} @ t={:.2f}".format(
                    arm, branch, b, len(proteins), f, t))

    with open(r / "summary.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["arm", "branch", "bin", "n_proteins", "fmax", "threshold"])
        w.writerows(rows)

    # ---- primary endpoint: bootstrap delta in lt30 / molecular_function ----
    branch, b = "molecular_function", "lt30"
    proteins = sorted(acc for acc in truth[branch] if bins.get(acc) == b)
    deltas = []
    if proteins:
        curves = {}
        for arm in ("A", "B"):
            curves[arm] = [
                protein_curve(
                    {t: s for t, s in arms[arm].get(acc, {}).items() if dag.branch(t) == branch},
                    truth[branch][acc],
                )
                for acc in proteins
            ]
        sizes = [len(truth[branch][acc]) for acc in proteins]
        rng = random.Random(cfg["seed"])
        for _ in range(N_BOOT):
            idx = [rng.randrange(len(proteins)) for _ in proteins]
            fa, _ = fmax_from_curves([curves["A"][i] for i in idx], [sizes[i] for i in idx])
            fb, _ = fmax_from_curves([curves["B"][i] for i in idx], [sizes[i] for i in idx])
            deltas.append(fb - fa)
        deltas.sort()
        lo = deltas[int(0.025 * N_BOOT)]
        hi = deltas[int(0.975 * N_BOOT) - 1]
        point = fmax_table[("B", branch, b)] - fmax_table[("A", branch, b)]
        verdict = "SUPPORTED" if lo > 0 else ("FALSIFIED" if hi < 0 else "INCONCLUSIVE")
    else:
        point, lo, hi, verdict = 0.0, 0.0, 0.0, "NO DATA IN BIN"

    with open(r / "report.md", "w", encoding="utf-8", newline="\n") as f:
        f.write("# v0.1 twilight-zone experiment - report\n\n")
        f.write("**Question:** does structure (AlphaFold + Foldseek) improve GO-term ")
        f.write("prediction where sequence similarity fails (<30% identity)?\n\n")
        f.write("**Primary endpoint** (molecular_function, lt30 bin, n={}):\n\n".format(len(proteins)))
        f.write("- delta Fmax (Arm B - Arm A): **{:+.4f}**  (95% bootstrap CI [{:+.4f}, {:+.4f}])\n".format(point, lo, hi))
        f.write("- Verdict: **{}**\n\n".format(verdict))
        f.write("## Fmax by arm / branch / identity bin\n\n")
        f.write("| branch | bin | n | Arm A (seq) | Arm B (+struct) | delta |\n")
        f.write("|---|---|---:|---:|---:|---:|\n")
        for branch2 in BRANCHES:
            for b2 in BIN_ORDER:
                n = len([a for a in truth[branch2] if bins.get(a) == b2])
                fa = fmax_table.get(("A", branch2, b2), 0.0)
                fb = fmax_table.get(("B", branch2, b2), 0.0)
                f.write("| {} | {} | {} | {:.4f} | {:.4f} | {:+.4f} |\n".format(
                    branch2, b2, n, fa, fb, fb - fa))
        f.write("\nEvidence type of every prediction here: COMPUTATIONAL. ")
        f.write("Ground truth: EXPERIMENTAL/CURATED GO annotations only.\n")

    print()
    print("PRIMARY ENDPOINT  delta Fmax (MF, lt30): {:+.4f}  CI [{:+.4f}, {:+.4f}]  -> {}".format(
        point, lo, hi, verdict))
    print("Wrote results/summary.csv and results/report.md")


if __name__ == "__main__":
    main()
