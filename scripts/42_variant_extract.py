"""Extract validated ClinVar missense variants for ALL human proteins in our
structure set (annotated + dark), not just dark ones.

Variants are validated by reference-residue match against our sequences
(rejects wrong isoforms). Output feeds the variant-at-interface analysis.

Outputs: results/explore/variants_human.tsv
         (accession, variant, significance, dark)
"""
import gzip
import re
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import data_dir, http_session, load_config, read_fasta
from pis.go import parse_gaf_full

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
       "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
       "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
       "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
P_RE = re.compile(r"\(p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})\)")
PATHO = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
BENIGN = {"Benign", "Likely benign", "Benign/Likely benign"}


def symbol_map(session, accs):
    m = defaultdict(list)
    accs = sorted(accs)
    for i in range(0, len(accs), 500):
        chunk = accs[i:i + 500]
        for attempt in range(3):
            try:
                r = session.post("https://mygene.info/v3/query",
                                 data={"q": ",".join(chunk), "scopes": "uniprot",
                                       "fields": "symbol", "species": "human"},
                                 timeout=120)
                r.raise_for_status()
                for row in r.json():
                    if row.get("symbol"):
                        m[row["symbol"]].append(row["query"])
                break
            except Exception:
                time.sleep(5 * (attempt + 1))
    return m


def main():
    cfg = load_config()
    d = data_dir(cfg)
    out = Path("results/explore/variants_human.tsv")
    session = http_session()

    kb = set((d / "final_accessions.txt").read_text().split())
    dark = set((d / "dark" / "accessions.txt").read_text().split())
    human = set()
    for acc, _t, _a, _e in parse_gaf_full(d / "gaf" / "goa_human.gaf.gz"):
        if acc in kb or acc in dark:
            human.add(acc)
    seqs = read_fasta(d / "final.fasta")
    seqs.update(read_fasta(d / "dark" / "dark.fasta"))
    print("human proteins in set: {}".format(len(human)))
    sym2acc = symbol_map(session, human)
    print("symbols mapped: {}".format(len(sym2acc)))

    rows = []
    with gzip.open(d / "explore" / "variant_summary.txt.gz", "rt",
                   encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in f:
            c = line.rstrip("\n").split("\t")
            if c[idx["Assembly"]] != "GRCh38":
                continue
            sig = c[idx["ClinicalSignificance"]].split(";")[0].strip()
            grp = "pathogenic" if sig in PATHO else ("benign" if sig in BENIGN else None)
            if not grp:
                continue
            gene = c[idx["GeneSymbol"]]
            if gene not in sym2acc:
                continue
            m = P_RE.search(c[idx["Name"]])
            if not m or m.group(1) not in AA3 or m.group(3) not in AA3:
                continue
            ref, pos, alt = AA3[m.group(1)], int(m.group(2)), AA3[m.group(3)]
            if ref == alt:
                continue
            for acc in sym2acc[gene]:
                seq = seqs.get(acc, "")
                if pos <= len(seq) and seq[pos - 1] == ref:
                    rows.append((acc, "{}{}{}".format(ref, pos, alt), grp,
                                 "1" if acc in dark else "0"))
                    break

    seen = set()
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        f.write("accession\tvariant\tsignificance\tdark\n")
        for acc, var, grp, dk in rows:
            key = (acc, var, grp)
            if key in seen:
                continue
            seen.add(key)
            f.write("{}\t{}\t{}\t{}\n".format(acc, var, grp, dk))
    np = sum(1 for r in rows if r[2] == "pathogenic")
    nb = sum(1 for r in rows if r[2] == "benign")
    prot = len({r[0] for r in rows if r[2] == "pathogenic"})
    print("validated: {} pathogenic, {} benign sites; {} proteins with a "
          "pathogenic variant".format(np, nb, prot))
    print("Wrote", out)


if __name__ == "__main__":
    main()
