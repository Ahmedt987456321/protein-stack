"""Step 05 - neighbour searches for function transfer.

Sequence arm:   MMseqs2  test.fasta  vs train.fasta   (similarity = identity)
Structure arm:  Foldseek test PDBs   vs train PDBs    (similarity = TM-score)

Outputs:
  data/seq_hits.tsv     (query, target, sim)   sim = identity fraction
  data/struct_hits.tsv  (query, target, sim)   sim = alignment TM-score
"""
import shutil
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import data_dir, load_config, norm_hit_id, run_tool


def normalize(raw: Path, out: Path, percent_ok: bool):
    n = 0
    with open(raw, encoding="utf-8") as fin, open(out, "w", encoding="utf-8", newline="\n") as fout:
        fout.write("query\ttarget\tsim\n")
        for line in fin:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 3:
                continue
            q, t = norm_hit_id(cols[0]), norm_hit_id(cols[1])
            if q == t:
                continue
            sim = float(cols[2])
            if percent_ok and sim > 1.0:
                sim /= 100.0
            fout.write("{}\t{}\t{:.4f}\n".format(q, t, sim))
            n += 1
    return n


def main():
    cfg = load_config()
    d = data_dir(cfg)
    ev = str(cfg["search"]["evalue"])
    ms = str(cfg["search"]["max_seqs"])
    th = str(cfg["tools"]["threads"])

    seq_raw = d / "seq_hits_raw.tsv"
    tmp1 = d / "tmp_seq"
    run_tool(cfg, "mmseqs", [
        "easy-search",
        d / "test.fasta",
        d / "train.fasta",
        seq_raw,
        tmp1,
        "--format-output", "query,target,pident,bits,evalue",
        "-e", ev, "-s", "7.5", "--max-seqs", ms, "--threads", th,
    ])
    n = normalize(seq_raw, d / "seq_hits.tsv", percent_ok=True)
    print("Sequence hits: {} -> data/seq_hits.tsv".format(n))
    shutil.rmtree(tmp1, ignore_errors=True)

    struct_raw = d / "struct_hits_raw.tsv"
    tmp2 = d / "tmp_struct"
    run_tool(cfg, "foldseek", [
        "easy-search",
        d / "struct_test",
        d / "struct_train",
        struct_raw,
        tmp2,
        "--format-output", "query,target,alntmscore,bits,evalue",
        "-e", ev, "--max-seqs", ms, "--threads", th,
    ])
    n = normalize(struct_raw, d / "struct_hits.tsv", percent_ok=False)
    print("Structure hits: {} -> data/struct_hits.tsv".format(n))
    shutil.rmtree(tmp2, ignore_errors=True)


if __name__ == "__main__":
    main()
