"""Step 20 (v1.0) - the prospective loop: re-check live GOA against the
frozen 2026 dark hypotheses.

Downloads today's GAFs (into data/gaf_live/, separate from the snapshot the
pipeline was built on), finds dark-set proteins that have gained EXPERIMENTAL
annotations since our snapshot, and scores the frozen hypotheses in
results/dark_hypotheses.tsv against them. Outcomes are appended to
results/hypothesis_outcomes.tsv (deduplicated) so step 13 folds them into the
graph on the next rebuild.

Run this any time - monthly is plenty (GOA releases are roughly quarterly).
Finding nothing new is the expected result at first.
"""
import datetime
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import data_dir, download, http_session, load_config, results_dir
from pis.go import GoDag, parse_gaf


def main():
    cfg = load_config()
    d = data_dir(cfg)
    live = d / "gaf_live"
    live.mkdir(exist_ok=True)
    r = results_dir(cfg)
    session = http_session()
    dag = GoDag(d / "go-basic.obo")
    today = datetime.date.today().isoformat()

    dark = set((d / "dark" / "accessions.txt").read_text(encoding="utf-8").split())

    # live experimental annotations for dark proteins
    live_exp = defaultdict(set)
    for sp in cfg["species"]:
        dest = live / "goa_{}_{}.gaf.gz".format(sp["name"], today)
        print("Fetching live GAF: {} ...".format(sp["name"]))
        download(session, sp["gaf_url"], dest)
        for acc, term, _aspect in parse_gaf(dest, cfg["evidence_codes"]):
            if acc in dark:
                live_exp[acc].add(term)

    newly = sorted(live_exp)
    print("Dark proteins with experimental annotations in live GOA: {}".format(len(newly)))
    if not newly:
        print("No validation events yet - expected shortly after the snapshot. "
              "Re-run after the next GOA release.")
        return

    # frozen hypotheses for those proteins
    hyps = defaultdict(list)
    with open(r / "dark_hypotheses.tsv", encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, branch, term, score, _tier, _streams, _so = line.rstrip("\n").split("\t")
            if acc in live_exp:
                hyps[acc].append((term, branch, float(score)))

    out_path = r / "hypothesis_outcomes.tsv"
    seen = set()
    if out_path.exists():
        with open(out_path, encoding="utf-8") as f:
            next(f)
            seen = {tuple(line.split("\t")[:2]) for line in f}
    else:
        with open(out_path, "w", encoding="utf-8", newline="\n") as f:
            f.write("accession\tgo_term\tbranch\tconfidence\tarm\toutcome\tsource\n")

    source = "GOA-prospective-{}".format(today)
    n_sup = n_unc = 0
    with open(out_path, "a", encoding="utf-8", newline="\n") as f:
        for acc in newly:
            truth = dag.propagate(live_exp[acc])
            for term, branch, score in hyps.get(acc, []):
                if (acc, term) in seen:
                    continue
                outcome = "supported" if term in truth else "unconfirmed"
                if outcome == "supported":
                    n_sup += 1
                else:
                    n_unc += 1
                f.write("{}\t{}\t{}\t{:.4f}\tFROZEN\t{}\t{}\n".format(
                    acc, term, branch, score, outcome, source))

    print("Prospective outcomes appended: supported={} unconfirmed={} "
          "(re-run step 13 to fold into kg.db)".format(n_sup, n_unc))


if __name__ == "__main__":
    main()
