"""Step 03 - download AlphaFold models, apply the mean-pLDDT quality gate,
fix the final protein set.

Uses the AlphaFold DB API per accession (robust to model-version bumps).
Files are saved as data/structures/<ACC>.pdb so Foldseek hit names map
straight back to accessions.

Outputs:
  data/structures/<ACC>.pdb   (only models passing the gate)
  data/structures.csv         (accession, mean_plddt, kept)
  data/final_accessions.txt   final protein set (<= max_proteins)
  data/final.fasta
"""
import csv
import random
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tqdm import tqdm

from pis.common import data_dir, http_session, load_config, read_fasta, write_fasta
from pis.fetch import fetch_alphafold as fetch_one


def main():
    cfg = load_config()
    d = data_dir(cfg)
    out_dir = d / "structures"
    out_dir.mkdir(exist_ok=True)
    session = http_session()

    seqs = read_fasta(d / "sequences.fasta")
    accs = sorted(seqs)
    min_plddt = cfg["structures"]["min_plddt"]
    print("Fetching AlphaFold models for {} proteins (gate: mean pLDDT >= {}) ...".format(len(accs), min_plddt))

    rows = []
    with ThreadPoolExecutor(max_workers=cfg["structures"]["workers"]) as pool:
        futures = [pool.submit(fetch_one, session, a, out_dir, min_plddt) for a in accs]
        for fut in tqdm(as_completed(futures), total=len(futures), unit="protein"):
            rows.append(fut.result())

    rows.sort()
    with open(d / "structures.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["accession", "mean_plddt", "kept"])
        for acc, plddt, kept in rows:
            w.writerow([acc, "{:.2f}".format(plddt), int(kept)])

    kept = [acc for acc, _, k in rows if k]
    print("Models kept after pLDDT gate: {} / {}".format(len(kept), len(accs)))

    cap = cfg["sampling"]["max_proteins"]
    if len(kept) > cap:
        rng = random.Random(cfg["seed"])
        rng.shuffle(kept)
        kept = sorted(kept[:cap])
        print("Downsampled to max_proteins = {}".format(cap))
    else:
        kept = sorted(kept)

    with open(d / "final_accessions.txt", "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(kept) + "\n")
    write_fasta(d / "final.fasta", {a: seqs[a] for a in kept})
    print("Final protein set: {} -> data/final_accessions.txt, data/final.fasta".format(len(kept)))


if __name__ == "__main__":
    main()
