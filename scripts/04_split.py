"""Step 04 - train/test split + identity bins.

Random protein-level split (seeded), then each test protein is binned by its
maximum sequence identity to the training set (MMseqs2 search, permissive
e-value so weak homology is still detected). No hit at all -> lt30 bin.

Outputs:
  data/train.fasta, data/test.fasta
  data/struct_train/, data/struct_test/   (hardlinked PDBs for Foldseek)
  data/split.csv                          (accession, split, bin)
"""
import csv
import random
import shutil
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import (
    data_dir,
    link_or_copy,
    load_config,
    norm_hit_id,
    read_fasta,
    run_tool,
    write_fasta,
)

BIN_LABELS = ["lt30", "30to50", "50to80", "ge80"]


def bin_of(identity, edges):
    """identity as fraction 0..1; edges e.g. [0.3, 0.5, 0.8]."""
    for i, edge in enumerate(edges):
        if identity < edge:
            return BIN_LABELS[i]
    return BIN_LABELS[len(edges)]


def main():
    cfg = load_config()
    d = data_dir(cfg)
    seqs = read_fasta(d / "final.fasta")
    accs = sorted(seqs)

    rng = random.Random(cfg["seed"])
    rng.shuffle(accs)
    n_test = round(cfg["split"]["test_frac"] * len(accs))
    test = sorted(accs[:n_test])
    train = sorted(accs[n_test:])
    print("Split: {} train / {} test".format(len(train), len(test)))

    write_fasta(d / "train.fasta", {a: seqs[a] for a in train})
    write_fasta(d / "test.fasta", {a: seqs[a] for a in test})

    for split_name, members in (("struct_train", train), ("struct_test", test)):
        out = d / split_name
        if out.exists():
            shutil.rmtree(out)
        out.mkdir()
        missing = 0
        for acc in members:
            src = d / "structures" / (acc + ".pdb")
            if src.exists():
                link_or_copy(src, out / (acc + ".pdb"))
            else:
                missing += 1
        print("{}: {} structures linked ({} missing)".format(split_name, len(members) - missing, missing))

    # max identity of each test protein to the training set
    tmp = d / "tmp_bins"
    hits = d / "bin_hits.tsv"
    run_tool(cfg, "mmseqs", [
        "easy-search",
        d / "test.fasta",
        d / "train.fasta",
        hits,
        tmp,
        "--format-output", "query,target,pident,evalue",
        "-e", "10",
        "-s", "7.5",
        "--max-seqs", "500",
        "--threads", str(cfg["tools"]["threads"]),
    ])

    max_id = defaultdict(float)
    with open(hits, encoding="utf-8") as f:
        for line in f:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            q = norm_hit_id(cols[0])
            ident = float(cols[2])
            if ident > 1.0:  # mmseqs reports percent
                ident /= 100.0
            max_id[q] = max(max_id[q], ident)

    edges = cfg["split"]["bin_edges"]
    counts = defaultdict(int)
    with open(d / "split.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["accession", "split", "bin"])
        for acc in train:
            w.writerow([acc, "train", ""])
        for acc in test:
            b = bin_of(max_id.get(acc, 0.0), edges)
            counts[b] += 1
            w.writerow([acc, "test", b])

    shutil.rmtree(tmp, ignore_errors=True)
    print("Test bins: " + ", ".join("{}={}".format(b, counts[b]) for b in BIN_LABELS))
    print("Wrote data/split.csv")


if __name__ == "__main__":
    main()
