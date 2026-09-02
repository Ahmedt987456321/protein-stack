"""Step 10 (v0.2) - evaluate all arms and check the phase-0.2 gates.

Arms:
  A  sequence-only transfer            (v0.1 baseline)
  B  naive max-fusion seq+structure    (v0.1 treatment)
  S  structure-only transfer
  D  domain-only (InterPro)
  C  fitted identity-aware fusion      (v0.2)

Gates:
  G1  C > A in the lt30 bin, molecular_function (bootstrap 95% CI > 0)
  G2  no statistically significant dilution: in every branch x bin cell where
      the C-A point estimate is negative, the bootstrap 95% CI must include 0
      (small cells have Fmax noise of several points; a point-estimate check
      would fail on noise alone)
  G3  C beats every single stream (A, S, D) on macro-average Fmax over cells

Outputs:
  results/summary_v2.csv
  results/report_v2.md
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

ARMS = ["A", "B", "S", "D", "C"]
ARM_FILES = {
    "A": "pred_armA.tsv", "B": "pred_armB.tsv", "S": "pred_armS.tsv",
    "D": "pred_armD.tsv", "C": "pred_armC.tsv",
}
ARM_NAMES = {
    "A": "sequence only", "B": "naive max-fusion (v0.1)",
    "S": "structure only", "D": "domains only", "C": "fitted fusion (v0.2)",
}
BRANCHES = ["molecular_function", "biological_process", "cellular_component"]
BIN_ORDER = ["lt30", "30to50", "50to80", "ge80"]
N_BOOT = 1000


def load_preds(path):
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

    bins = {}
    with open(d / "split.csv", encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, split, b = line.rstrip("\n").split(",")
            if split == "test":
                bins[acc] = b

    raw_ann = defaultdict(set)
    with open(d / "annotations.tsv", encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, term, _aspect = line.rstrip("\n").split("\t")
            if acc in bins:
                raw_ann[acc].add(term)

    truth = {b: {} for b in BRANCHES}
    for acc, terms in raw_ann.items():
        prop = dag.propagate(terms)
        for b in BRANCHES:
            sub = {t for t in prop if dag.branch(t) == b}
            if sub:
                truth[b][acc] = sub

    arms = {a: load_preds(r / ARM_FILES[a]) for a in ARMS}

    # per-branch, per-arm branch-filtered predictions (built once)
    branch_preds = {}
    for arm in ARMS:
        for branch in BRANCHES:
            branch_preds[(arm, branch)] = {
                acc: {t: s for t, s in arms[arm].get(acc, {}).items() if dag.branch(t) == branch}
                for acc in truth[branch]
            }

    table = {}   # (arm, branch, bin) -> (fmax, n)
    rows = []
    for branch in BRANCHES:
        for b in BIN_ORDER:
            proteins = sorted(acc for acc in truth[branch] if bins.get(acc) == b)
            for arm in ARMS:
                f, t = fmax(proteins, branch_preds[(arm, branch)], truth[branch])
                table[(arm, branch, b)] = (f, len(proteins))
                rows.append([arm, branch, b, len(proteins), "{:.4f}".format(f), "{:.2f}".format(t)])

    with open(r / "summary_v2.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["arm", "branch", "bin", "n_proteins", "fmax", "threshold"])
        w.writerows(rows)

    cells = [(branch, b) for branch in BRANCHES for b in BIN_ORDER
             if table[("A", branch, b)][1] > 0]
    macro = {arm: sum(table[(arm, br, b)][0] for br, b in cells) / len(cells) for arm in ARMS}

    def bootstrap_ca(branch, b, rng):
        """Bootstrap 95% CI for Fmax(C) - Fmax(A) in one cell."""
        proteins = sorted(acc for acc in truth[branch] if bins.get(acc) == b)
        curves = {
            arm: [protein_curve(branch_preds[(arm, branch)].get(acc, {}), truth[branch][acc])
                  for acc in proteins]
            for arm in ("A", "C")
        }
        sizes = [len(truth[branch][acc]) for acc in proteins]
        deltas = []
        for _ in range(N_BOOT):
            idx = [rng.randrange(len(proteins)) for _ in proteins]
            fa, _ = fmax_from_curves([curves["A"][i] for i in idx], [sizes[i] for i in idx])
            fc, _ = fmax_from_curves([curves["C"][i] for i in idx], [sizes[i] for i in idx])
            deltas.append(fc - fa)
        deltas.sort()
        return deltas[int(0.025 * N_BOOT)], deltas[int(0.975 * N_BOOT) - 1], len(proteins)

    rng = random.Random(cfg["seed"])

    # ---- G1: bootstrap C-A in (MF, lt30) -----------------------------------
    branch, b = "molecular_function", "lt30"
    lo, hi, n_primary = bootstrap_ca(branch, b, rng)
    point = table[("C", branch, b)][0] - table[("A", branch, b)][0]
    g1 = lo > 0

    # ---- G2: no statistically significant dilution --------------------------
    # bootstrap only the cells with a negative point estimate
    dilution_checks = []  # (branch, bin, point_delta, lo, hi, significant)
    for br, bb in cells:
        delta_pt = table[("C", br, bb)][0] - table[("A", br, bb)][0]
        if delta_pt < 0:
            clo, chi, _n = bootstrap_ca(br, bb, rng)
            dilution_checks.append((br, bb, delta_pt, clo, chi, chi < 0))
    diluted = [(br, bb) for br, bb, _d, _l, _h, sig in dilution_checks if sig]
    g2 = not diluted

    # ---- G3: fusion beats every single stream on macro-average -------------
    g3 = all(macro["C"] > macro[a] for a in ("A", "S", "D"))

    verdict = "PASSED" if (g1 and g2 and g3) else "NOT PASSED"

    with open(r / "report_v2.md", "w", encoding="utf-8", newline="\n") as f:
        f.write("# v0.2 - domains + identity-aware fusion - report\n\n")
        f.write("**Arms:** " + "; ".join("{} = {}".format(a, ARM_NAMES[a]) for a in ARMS) + "\n\n")
        f.write("## Gates\n\n")
        f.write("- **G1** (fusion beats sequence in the twilight zone, MF/lt30, n={}): "
                "delta Fmax {:+.4f}, 95% CI [{:+.4f}, {:+.4f}] -> {}\n".format(
                    n_primary, point, lo, hi, "PASS" if g1 else "FAIL"))
        f.write("- **G2** (no statistically significant dilution, per-cell bootstrap 95% CI): {}\n".format(
            "PASS" if g2 else "FAIL"))
        if dilution_checks:
            for br, bb, delta_pt, clo, chi, sig in dilution_checks:
                f.write("  - {}/{}: point {:+.4f}, CI [{:+.4f}, {:+.4f}] -> {}\n".format(
                    br, bb, delta_pt, clo, chi,
                    "SIGNIFICANT dilution" if sig else "not significant (CI includes 0)"))
        else:
            f.write("  - no cell had a negative point estimate\n")
        f.write("- **G3** (fusion beats every single stream, macro-Fmax): {}\n".format(
            "PASS" if g3 else "FAIL"))
        f.write("  - macro-Fmax: " + ", ".join(
            "{} {:.4f}".format(a, macro[a]) for a in ARMS) + "\n")
        f.write("\n**Phase 0.2 gate: {}**\n\n".format(verdict))
        f.write("## Fmax by branch / bin\n\n")
        f.write("| branch | bin | n | A seq | B naive | S struct | D domain | C fused | C-A |\n")
        f.write("|---|---|---:|---:|---:|---:|---:|---:|---:|\n")
        for br in BRANCHES:
            for bb in BIN_ORDER:
                fa = table[("A", br, bb)][0]
                fc = table[("C", br, bb)][0]
                f.write("| {} | {} | {} |".format(br, bb, table[("A", br, bb)][1]))
                for arm in ARMS:
                    f.write(" {:.4f} |".format(table[(arm, br, bb)][0]))
                f.write(" {:+.4f} |\n".format(fc - fa))
        f.write("\nEvidence types: seq/structure transfer and learned domain "
                "associations are COMPUTATIONAL; interpro2go mappings are CURATED; "
                "ground truth is EXPERIMENTAL GO only.\n")

    for br in BRANCHES:
        for bb in BIN_ORDER:
            print("{:22s} {:7s} n={:5d}  ".format(br, bb, table[("A", br, bb)][1]) +
                  "  ".join("{}={:.4f}".format(a, table[(a, br, bb)][0]) for a in ARMS))
    print()
    print("macro-Fmax: " + ", ".join("{}={:.4f}".format(a, macro[a]) for a in ARMS))
    print("G1 (MF/lt30 C-A): {:+.4f} CI [{:+.4f}, {:+.4f}] -> {}".format(
        point, lo, hi, "PASS" if g1 else "FAIL"))
    for br, bb, delta_pt, clo, chi, sig in dilution_checks:
        print("G2 check {}/{}: point {:+.4f} CI [{:+.4f}, {:+.4f}] -> {}".format(
            br, bb, delta_pt, clo, chi, "SIGNIFICANT" if sig else "not significant"))
    print("G2 (no significant dilution): {}".format("PASS" if g2 else "FAIL - " + str(diluted)))
    print("G3 (beats single streams): {}".format("PASS" if g3 else "FAIL"))
    print("PHASE 0.2 GATE: " + verdict)
    print("Wrote results/summary_v2.csv and results/report_v2.md")


if __name__ == "__main__":
    main()
