"""Step 14 (v0.4) - STRING interaction edges for our proteins.

Downloads per-species STRING links + aliases (bulk files), maps STRING ids to
UniProt accessions, and keeps high-confidence edges where BOTH ends are in our
protein universe (annotated KB + dark set).

Outputs:
  data/string/<taxid>.links.txt.gz, <taxid>.aliases.txt.gz  (cached downloads)
  data/string_edges.tsv   (acc_a, acc_b, combined_score 0-1000), a < b
"""
import gzip
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import data_dir, download, http_session, load_config

LINKS_URL = "https://stringdb-downloads.org/download/protein.links.v{v}/{tax}.protein.links.v{v}.txt.gz"
ALIAS_URL = "https://stringdb-downloads.org/download/protein.aliases.v{v}/{tax}.protein.aliases.v{v}.txt.gz"


def main():
    cfg = load_config()
    d = data_dir(cfg)
    sd = d / "string"
    sd.mkdir(exist_ok=True)
    session = http_session()
    v = cfg["string"]["version"]
    score_min = cfg["string"]["score_min"]

    ours = set((d / "final_accessions.txt").read_text(encoding="utf-8").split())
    ours |= set((d / "dark" / "accessions.txt").read_text(encoding="utf-8").split())
    print("Protein universe: {}".format(len(ours)))

    edges = {}
    for sp in cfg["species"]:
        tax = sp["taxid"]
        try:
            links = download(session, LINKS_URL.format(v=v, tax=tax), sd / "{}.links.txt.gz".format(tax))
            alias = download(session, ALIAS_URL.format(v=v, tax=tax), sd / "{}.aliases.txt.gz".format(tax))
        except Exception as e:
            print("  warn: STRING unavailable for {} (taxid {}): {}; skipping".format(
                sp["name"], tax, e))
            continue

        s2u = {}
        with gzip.open(alias, "rt", encoding="utf-8", errors="replace") as f:
            for line in f:
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 3 or "UniProt_AC" not in cols[2]:
                    continue
                acc = cols[1].split("-")[0]
                if acc in ours:
                    s2u[cols[0]] = acc
        print("{}: {} STRING ids map to our proteins".format(sp["name"], len(s2u)))

        n = 0
        with gzip.open(links, "rt", encoding="utf-8", errors="replace") as f:
            next(f)  # header
            for line in f:
                cols = line.rstrip("\n").split()
                if len(cols) < 3:
                    continue
                a = s2u.get(cols[0])
                b = s2u.get(cols[1])
                if a is None or b is None or a == b:
                    continue
                score = int(cols[2])
                if score < score_min:
                    continue
                key = (a, b) if a < b else (b, a)
                if score > edges.get(key, 0):
                    edges[key] = score
                    n += 1
        print("{}: {} qualifying edge rows".format(sp["name"], n))

    with open(d / "string_edges.tsv", "w", encoding="utf-8", newline="\n") as f:
        f.write("acc_a\tacc_b\tcombined_score\n")
        for (a, b), score in sorted(edges.items()):
            f.write("{}\t{}\t{}\n".format(a, b, score))
    print("STRING edges (score >= {}): {} -> data/string_edges.tsv".format(score_min, len(edges)))


if __name__ == "__main__":
    main()
