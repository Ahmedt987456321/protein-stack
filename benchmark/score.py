"""timesplit-go scorer. Self-contained; requires only go-basic.obo.

Usage:
  python score.py --predictions preds.tsv --truth truth_2023.tsv \
      --obo go-basic.obo [--permutations 1000]

Predictions file: tab-separated, one row per protein and branch:
  accession <TAB> branch <TAB> go_term
where branch is molecular_function, biological_process, or
cellular_component and go_term is the submitter's single most specific
prediction for that protein and branch. Proteins or branches may be
omitted; only submitted rows are scored (coverage is reported).

Scoring: a prediction is a hit if the predicted term is in the truth set
for that protein and branch (truth is already propagated over is_a and
part_of, so predicting a correct ancestor scores as a hit). Reported per
branch: n scored, top-1 precision, and a permutation p-value obtained by
shuffling the submitted terms across proteins.
"""
import argparse
import random
from collections import defaultdict

BRANCHES = ("molecular_function", "biological_process", "cellular_component")


def load_obo(path):
    parents, alt, obsolete = {}, {}, set()
    term = None
    in_term = False
    with open(path, encoding="utf-8") as f:
        for raw in f:
            line = raw.strip()
            if line == "[Term]":
                in_term, term = True, None
            elif line.startswith("["):
                in_term = False
            elif in_term and line.startswith("id: "):
                term = line[4:]
                parents.setdefault(term, set())
            elif in_term and term:
                if line.startswith("alt_id: "):
                    alt[line[8:]] = term
                elif line.startswith("is_a: "):
                    parents[term].add(line[6:].split(" ")[0])
                elif line.startswith("relationship: part_of "):
                    parents[term].add(line[22:].split(" ")[0])
                elif line == "is_obsolete: true":
                    obsolete.add(term)
    return parents, alt, obsolete


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--predictions", required=True)
    ap.add_argument("--truth", required=True)
    ap.add_argument("--obo", required=True)
    ap.add_argument("--ic", help="term IC table (term<TAB>ic); enables the"
                    " primary information-gain metric")
    ap.add_argument("--permutations", type=int, default=1000)
    ap.add_argument("--seed", type=int, default=17)
    args = ap.parse_args()

    parents, alt, obsolete = load_obo(args.obo)

    truth = defaultdict(lambda: defaultdict(set))
    with open(args.truth, encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, branch, term = line.rstrip("\n").split("\t")
            truth[branch][acc].add(term)

    preds = defaultdict(dict)
    with open(args.predictions, encoding="utf-8") as f:
        first = f.readline()
        if not first.lower().startswith("accession"):
            f.seek(0)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            acc, branch, term = parts[0], parts[1], alt.get(parts[2], parts[2])
            if branch in BRANCHES and term not in obsolete:
                preds[branch][acc] = term

    ic = {}
    if args.ic:
        with open(args.ic, encoding="utf-8") as f:
            next(f)
            for line in f:
                t, v = line.rstrip("\n").split("\t")
                ic[t] = float(v)
    ic_default = max(ic.values()) if ic else 0.0

    rng = random.Random(args.seed)
    print("branch\tn_truth\tn_scored\ttop1_precision\tperm_p\tmean_info_gain_bits")
    for branch in BRANCHES:
        cohort = truth[branch]
        scored = sorted(a for a in cohort if a in preds[branch])
        if not scored:
            print("{}\t{}\t0\t-\t-".format(branch, len(cohort)))
            continue
        hits = sum(1 for a in scored if preds[branch][a] in cohort[a])
        obs = hits / len(scored)
        terms = [preds[branch][a] for a in scored]
        ge = 0
        for _ in range(args.permutations):
            rng.shuffle(terms)
            h = sum(1 for a, t in zip(scored, terms) if t in cohort[a])
            if h / len(scored) >= obs:
                ge += 1
        p = (ge + 1) / (args.permutations + 1)
        if ic:
            gain = sum(ic.get(preds[branch][a], ic_default)
                       for a in scored if preds[branch][a] in cohort[a])
            mig = "{:.3f}".format(gain / len(scored))
        else:
            mig = "-"
        print("{}\t{}\t{}\t{:.4f}\t{:.4f}\t{}".format(
            branch, len(cohort), len(scored), obs, p, mig))


if __name__ == "__main__":
    main()
