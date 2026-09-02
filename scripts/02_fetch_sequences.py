"""Step 02 - fetch canonical sequences from UniProt for the sampled accessions,
apply the length filter.

Outputs:
  data/sequences.fasta   (>ACC headers, length-filtered)
"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tqdm import tqdm

from pis.common import data_dir, http_session, load_config, write_fasta
from pis.fetch import UNIPROT_BATCH as BATCH
from pis.fetch import fetch_uniprot_batch as fetch_batch


def main():
    cfg = load_config()
    d = data_dir(cfg)
    session = http_session()

    accs = (d / "accessions.txt").read_text(encoding="utf-8").split()
    print("Fetching sequences for {} accessions ...".format(len(accs)))

    seqs = {}
    for i in tqdm(range(0, len(accs), BATCH), unit="batch"):
        seqs.update(fetch_batch(session, accs[i : i + BATCH]))

    lo, hi = cfg["sampling"]["min_len"], cfg["sampling"]["max_len"]
    kept = {a: s for a, s in seqs.items() if lo <= len(s) <= hi}
    print("Fetched {}; kept {} after length filter [{}, {}]".format(len(seqs), len(kept), lo, hi))

    write_fasta(d / "sequences.fasta", kept)
    print("Wrote data/sequences.fasta")


if __name__ == "__main__":
    main()
