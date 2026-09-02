"""Exploration B3 - ClinVar missense variants on our human structures.

Missense variants (GRCh38 rows) are mapped to our human proteins by gene
symbol, then validated: the variant's reference amino acid must match our
sequence at that position (this rejects wrong isoforms and stale symbols).
Each validated variant site is annotated with its pLDDT; for dark proteins,
fpocket is re-run on the affected subset to test whether pathogenic sites
fall inside druggable pockets.

Outputs: results/explore/variants.md, variants_dark.tsv
"""
import gzip
import re
import subprocess
import sys
import time
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import data_dir, http_session, load_config, read_fasta, to_wsl_path
from pis.go import parse_gaf_full

AA3 = {"Ala": "A", "Arg": "R", "Asn": "N", "Asp": "D", "Cys": "C",
       "Gln": "Q", "Glu": "E", "Gly": "G", "His": "H", "Ile": "I",
       "Leu": "L", "Lys": "K", "Met": "M", "Phe": "F", "Pro": "P",
       "Ser": "S", "Thr": "T", "Trp": "W", "Tyr": "Y", "Val": "V"}
P_RE = re.compile(r"\(p\.([A-Za-z]{3})(\d+)([A-Za-z]{3})\)")
PATHO = {"Pathogenic", "Likely pathogenic", "Pathogenic/Likely pathogenic"}
BENIGN = {"Benign", "Likely benign", "Benign/Likely benign"}


def symbol_map(session, accs):
    mapping = defaultdict(list)
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
                    sym = row.get("symbol")
                    if sym:
                        mapping[sym].append(row["query"])
                break
            except Exception:
                time.sleep(5 * (attempt + 1))
    return mapping


def residue_plddt(pdb_path):
    vals = {}
    with open(pdb_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    vals[int(line[22:26])] = float(line[60:66])
                except ValueError:
                    pass
    return vals


def pocket_residues_wsl(cfg, d, accs):
    """Re-run fpocket on the given dark accessions; return acc -> set of
    residue numbers in pockets with druggability >= 0.5."""
    workdir = d / "explore" / "pocket_res"
    workdir.mkdir(parents=True, exist_ok=True)
    listing = workdir / "accs.txt"
    with open(listing, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(accs) + "\n")
    script = workdir / "run.sh"
    script_text = (
        "set -e\nSRC={}\nOUT={}\nrm -rf /tmp/pr && mkdir -p /tmp/pr\n"
        "while read a; do cp \"$SRC/$a.pdb\" /tmp/pr/ 2>/dev/null || true; done < {}\n"
        "cd /tmp/pr\nls *.pdb | xargs -P {} -I XX fpocket -f XX > /dev/null 2>&1 || true\n"
        ": > residues.tsv\n: > goodpockets.tsv\n"
        "for info in *_out/*_info.txt; do\n"
        "  acc=$(basename \"$info\" _info.txt)\n"
        "  awk -v acc=\"$acc\" '/^Pocket/ {{pk=$2}} /Druggability Score/ "
        "{{d[pk]=$4+0}} END {{for (k in d) if (d[k]>=0.5) print acc\"\\t\"k}}' "
        "\"$info\" >> goodpockets.tsv\ndone\n"
        "while IFS=$'\\t' read acc pk; do\n"
        "  f=\"${{acc}}_out/pockets/pocket${{pk}}_atm.pdb\"\n"
        "  [ -f \"$f\" ] && awk -v acc=\"$acc\" '/^ATOM/ "
        "{{print acc\"\\t\"substr($0,23,4)+0}}' \"$f\" >> residues.tsv\n"
        "done < goodpockets.tsv\n"
        "sort -u residues.tsv > \"$OUT\"\nrm -rf /tmp/pr\n").format(
            to_wsl_path(d / "dark" / "structures"),
            to_wsl_path(workdir / "pocket_residues.tsv"),
            to_wsl_path(listing), cfg["tools"]["threads"])
    with open(script, "w", encoding="utf-8", newline="\n") as f:
        f.write(script_text)
    subprocess.run(["wsl", "-u", "root", "bash", to_wsl_path(script)], check=True)
    out = defaultdict(set)
    with open(workdir / "pocket_residues.tsv", encoding="utf-8") as f:
        for line in f:
            acc, res = line.rstrip("\n").split("\t")
            out[acc].add(int(res))
    return out


def main():
    cfg = load_config()
    d = data_dir(cfg)
    out_dir = Path("results/explore")
    out_dir.mkdir(parents=True, exist_ok=True)
    session = http_session()

    dark = set((d / "dark" / "accessions.txt").read_text(encoding="utf-8").split())
    kb = set((d / "final_accessions.txt").read_text(encoding="utf-8").split())
    human = set()
    for acc, _t, _a, _e in parse_gaf_full(d / "gaf" / "goa_human.gaf.gz"):
        if acc in dark or acc in kb:
            human.add(acc)
    seqs = read_fasta(d / "final.fasta")
    seqs.update(read_fasta(d / "dark" / "dark.fasta"))

    print("mapping {} human proteins to symbols ...".format(len(human)))
    sym2acc = symbol_map(session, human)
    print("symbols mapped: {}".format(len(sym2acc)))

    # ---- stream ClinVar ----------------------------------------------------
    counts = defaultdict(int)
    variants = []  # (acc, pos, ref, alt, significance_group, is_dark)
    with gzip.open(d / "explore" / "variant_summary.txt.gz", "rt",
                   encoding="utf-8", errors="replace") as f:
        header = f.readline().rstrip("\n").split("\t")
        idx = {c: i for i, c in enumerate(header)}
        for line in f:
            c = line.rstrip("\n").split("\t")
            if c[idx["Assembly"]] != "GRCh38":
                continue
            sig = c[idx["ClinicalSignificance"]].split(";")[0].strip()
            if sig in PATHO:
                group = "pathogenic"
            elif sig in BENIGN:
                group = "benign"
            else:
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
                    variants.append((acc, pos, ref, alt, group, acc in dark))
                    counts[(group, acc in dark)] += 1
                    break
    print("validated missense variants:", dict(counts))

    # ---- pLDDT at variant sites -------------------------------------------
    plddt_cache = {}
    def site_plddt(acc, pos):
        if acc not in plddt_cache:
            for sd in (d / "structures", d / "dark" / "structures"):
                p = sd / (acc + ".pdb")
                if p.exists():
                    plddt_cache[acc] = residue_plddt(p)
                    break
            else:
                plddt_cache[acc] = {}
        return plddt_cache[acc].get(pos)

    dark_var_accs = sorted({v[0] for v in variants if v[5]})
    print("dark proteins with validated variants: {}".format(len(dark_var_accs)))
    pockets = pocket_residues_wsl(cfg, d, dark_var_accs) if dark_var_accs else {}

    import statistics
    groups = defaultdict(list)
    dark_rows = []
    for acc, pos, ref, alt, grp, is_dark in variants:
        pl = site_plddt(acc, pos)
        if pl is None:
            continue
        groups[(grp, is_dark)].append(pl)
        if is_dark:
            in_pocket = pos in pockets.get(acc, set())
            dark_rows.append([acc, "{}{}{}".format(ref, pos, alt), grp,
                              "{:.1f}".format(pl), "yes" if in_pocket else "no"])

    with open(out_dir / "variants_dark.tsv", "w", encoding="utf-8", newline="\n") as f:
        f.write("accession\tvariant\tsignificance\tsite_pLDDT\tin_druggable_pocket\n")
        for r in sorted(dark_rows):
            f.write("\t".join(r) + "\n")

    def med(k):
        v = groups.get(k, [])
        return (len(v), statistics.median(v) if v else 0.0)
    n_pd, m_pd = med(("pathogenic", True))
    n_bd, m_bd = med(("benign", True))
    n_pa, m_pa = med(("pathogenic", False))
    n_ba, m_ba = med(("benign", False))
    dark_patho_pocket = sum(1 for r in dark_rows
                            if r[2] == "pathogenic" and r[4] == "yes")
    dark_patho = sum(1 for r in dark_rows if r[2] == "pathogenic")
    dark_ben_pocket = sum(1 for r in dark_rows
                          if r[2] == "benign" and r[4] == "yes")
    dark_ben = sum(1 for r in dark_rows if r[2] == "benign")

    with open(out_dir / "variants.md", "w", encoding="utf-8", newline="\n") as f:
        f.write("# ClinVar missense variants on the dataset structures\n\n")
        f.write("Variants validated by reference-residue match against our "
                "sequences (wrong isoforms rejected).\n\n")
        f.write("| group | set | n sites | median site pLDDT |\n|---|---|---|---|\n")
        f.write("| pathogenic | dark | {} | {:.1f} |\n".format(n_pd, m_pd))
        f.write("| benign | dark | {} | {:.1f} |\n".format(n_bd, m_bd))
        f.write("| pathogenic | annotated | {} | {:.1f} |\n".format(n_pa, m_pa))
        f.write("| benign | annotated | {} | {:.1f} |\n\n".format(n_ba, m_ba))
        f.write("Dark proteins carrying pathogenic missense variants: "
                "{} variants across {} proteins; {} of {} pathogenic sites "
                "({:.0f}%) fall inside a druggable pocket, versus {} of {} "
                "benign sites ({:.0f}%).\n\n".format(
                    dark_patho, len({r[0] for r in dark_rows if r[2] == "pathogenic"}),
                    dark_patho_pocket, dark_patho,
                    100 * dark_patho_pocket / dark_patho if dark_patho else 0,
                    dark_ben_pocket, dark_ben,
                    100 * dark_ben_pocket / dark_ben if dark_ben else 0))
        f.write("Full per-variant table: variants_dark.tsv.\n")
    print("dark pathogenic sites: {} ({} in pockets); annotated pathogenic: {}".format(
        dark_patho, dark_patho_pocket, n_pa))
    print("Wrote results/explore/variants.md")


if __name__ == "__main__":
    main()
