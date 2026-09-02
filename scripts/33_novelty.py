"""Exploration 2 - structural novelty screen.

Searches every dark-set AlphaFold model against the experimentally
determined PDB (Foldseek prebuilt database, WSL-local) and combines with
the existing dark-versus-KB search. A protein whose best TM-score to both
is below 0.5 has no known structural neighbour and is a candidate novel
fold (or a poorly modelled region; mean pLDDT is reported alongside).

Outputs: results/explore/novelty.md, novelty.tsv
"""
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import data_dir, load_config, norm_hit_id, run_tool
from pis.fetch import mean_ca_plddt
from pis.go import GoDag
from pis.streams import load_domains

PDB_DB = "/root/dbs/pdb"  # WSL-local Foldseek database


def main():
    cfg = load_config()
    d = data_dir(cfg)
    out_dir = Path("results/explore")
    out_dir.mkdir(parents=True, exist_ok=True)
    dag = GoDag(d / "go-basic.obo")

    hits_path = d / "explore" / "dark_vs_pdb.tsv"
    if not hits_path.exists():
        run_tool(cfg, "foldseek", [
            "easy-search", d / "dark" / "structures", PDB_DB,
            hits_path, d / "explore" / "tmp_pdb_search",
            "--format-output", "query,target,alntmscore,evalue",
            "-e", "10", "--max-seqs", "50",
            "--threads", str(cfg["tools"]["threads"]),
        ])

    best_pdb = defaultdict(float)
    with open(hits_path, encoding="utf-8") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 3:
                q = norm_hit_id(c[0])
                best_pdb[q] = max(best_pdb[q], float(c[2]))

    best_kb = defaultdict(float)
    with open(d / "dark" / "struct_hits_raw.tsv", encoding="utf-8") as f:
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) >= 3:
                q = norm_hit_id(c[0])
                best_kb[q] = max(best_kb[q], float(c[2]))

    dark = (d / "dark" / "accessions.txt").read_text(encoding="utf-8").split()
    dark_domains = load_domains(d / "dark" / "domains.tsv")
    tiers = {}
    with open(Path("results") / "dark_hypotheses.tsv", encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, _b, _t, score, tier, _s, _so = line.rstrip("\n").split("\t")
            if acc not in tiers or tier == "HIGH":
                tiers[acc] = tier

    rows = []
    for acc in dark:
        p = best_pdb.get(acc, 0.0)
        k = best_kb.get(acc, 0.0)
        best = max(p, k)
        plddt = mean_ca_plddt(d / "dark" / "structures" / (acc + ".pdb"))
        rows.append([acc, "{:.3f}".format(p), "{:.3f}".format(k),
                     "{:.3f}".format(best), "{:.1f}".format(plddt),
                     str(len(dark_domains.get(acc, ()))), tiers.get(acc, "-")])
    rows.sort(key=lambda r: float(r[3]))

    with open(out_dir / "novelty.tsv", "w", encoding="utf-8", newline="\n") as f:
        f.write("accession\tbest_TM_PDB\tbest_TM_KB\tbest_TM_any\tmean_pLDDT\t"
                "n_interpro\thypothesis_tier\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    cands = [r for r in rows if float(r[3]) < 0.5]
    strong = [r for r in cands if float(r[4]) >= 80 and r[5] == "0"]
    with open(out_dir / "novelty.md", "w", encoding="utf-8", newline="\n") as f:
        f.write("# Structural novelty screen (dark set vs PDB and KB)\n\n")
        f.write("{} of {} dark proteins have no structural neighbour with "
                "TM >= 0.5 in either the PDB or the annotated KB; {} of "
                "those are high-confidence models (mean pLDDT >= 80) with "
                "no InterPro entry - the strongest candidate novel "
                "folds.\n\n".format(len(cands), len(rows), len(strong)))
        f.write("Top 25 by lowest best-TM:\n\n")
        f.write("| accession | TM vs PDB | TM vs KB | mean pLDDT | InterPro "
                "entries | tier |\n|---|---|---|---|---|---|\n")
        for r in cands[:25]:
            f.write("| {} | {} | {} | {} | {} | {} |\n".format(
                r[0], r[1], r[2], r[4], r[5], r[6]))
        f.write("\nCaveat: a low best-TM can reflect modelling failure as "
                "well as fold novelty; the pLDDT column separates the two "
                "readings.\n")
    print("novelty candidates (best TM < 0.5): {} | strong (pLDDT>=80, "
          "no InterPro): {}".format(len(cands), len(strong)))
    print("Wrote results/explore/novelty.md")


if __name__ == "__main__":
    main()
