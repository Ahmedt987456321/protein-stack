"""Step 01 - download GO ontology + GOA GAF files, keep experimental annotations,
sample candidate proteins.

Outputs:
  data/go-basic.obo
  data/gaf/goa_<species>.gaf.gz
  data/annotations.tsv    (accession, go_term, aspect) for sampled proteins
  data/accessions.txt     sampled accession list (oversampled; later steps filter)
"""
import random
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import data_dir, download, http_session, load_config
from pis.go import parse_gaf


def main():
    cfg = load_config()
    d = data_dir(cfg)
    session = http_session()

    print("Downloading GO ontology ...")
    download(session, cfg["go_obo_url"], d / "go-basic.obo")

    ann = defaultdict(set)  # acc -> {(term, aspect)}
    for sp in cfg["species"]:
        gaf = d / "gaf" / "goa_{}.gaf.gz".format(sp["name"])
        print("Downloading GAF: {} ...".format(sp["name"]))
        download(session, sp["gaf_url"], gaf)
        n = 0
        for acc, term, aspect in parse_gaf(gaf, cfg["evidence_codes"]):
            ann[acc].add((term, aspect))
            n += 1
        print("  {}: {} experimental annotation rows".format(sp["name"], n))

    accs = sorted(ann)
    print("Proteins with >=1 experimental GO term: {}".format(len(accs)))

    rng = random.Random(cfg["seed"])
    rng.shuffle(accs)
    sample = accs[: cfg["sampling"]["oversample"]]

    with open(d / "accessions.txt", "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(sorted(sample)) + "\n")

    with open(d / "annotations.tsv", "w", encoding="utf-8", newline="\n") as f:
        f.write("accession\tgo_term\taspect\n")
        for acc in sorted(sample):
            for term, aspect in sorted(ann[acc]):
                f.write("{}\t{}\t{}\n".format(acc, term, aspect))

    print("Sampled {} proteins -> data/accessions.txt, data/annotations.tsv".format(len(sample)))


if __name__ == "__main__":
    main()
