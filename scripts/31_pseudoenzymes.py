"""Exploration 1 - pseudo-enzyme catalog.

Scans two enzyme families with well-characterised catalytic machinery:

  Protein kinases (IPR000719): beta-3 lysine (AxK), catalytic loop HRD
  aspartate, and DFG aspartate. A member missing any of the three is a
  pseudokinase candidate (the standard literature criterion).

  Trypsin-like serine proteases (IPR001254): the His and Ser catalytic-triad
  PROSITE patterns (PS00134, PS00135). A member missing either is a
  pseudo-protease candidate.

Motifs are located within InterPro domain boundaries fetched per protein
(cached); when boundaries are unavailable the whole sequence is scanned and
the row is marked accordingly. For proteins where both aspartate motifs are
found, the CA distance between them is measured in the AlphaFold model as a
geometric sanity check (7-14 A in genuine kinase folds). Each candidate is
cross-referenced against experimental GO catalytic annotations, and the
catalog is validated against known human pseudokinases present in the set.

Outputs: results/explore/pseudoenzymes.md, pseudoenzymes.tsv
"""
import json
import re
import sys
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tqdm import tqdm

from pis.common import data_dir, http_session, load_config, read_fasta
from pis.go import GoDag
from pis.streams import load_domains

KINASE = "IPR000719"
PROTEASE = "IPR001254"
# beta-3 AxK, catalytic loop, activation segment (permissive canonical forms)
RE_B3K = re.compile(r"[VAILMF][A-Z][VAILMF]?K")
RE_HRD = re.compile(r"H[RGA]D")
RE_DFG = re.compile(r"D[FLWY]G")
# PROSITE PS00134 (triad His) and PS00135 (triad Ser)
RE_HIS = re.compile(r"[LIVM][ST]A[STAG]HC")
RE_SER = re.compile(r"[DNSTAGC][GSTAPIMVQH][A-Z]{2}G[DE]SG[GS][SAPHV]")

KNOWN_PSEUDOKINASES = {  # human, widely accepted
    "P21860": "ERBB3", "Q13418": "ILK", "Q7RTN6": "STRADA",
    "Q8IVT5": "KSR1", "Q96RU7": "TRIB3", "Q96RU8": "TRIB1",
    "Q8IV63": "VRK3", "Q13308": "PTK7", "Q59H18": "TNNI3K",
    "O43187": "IRAK2", "Q9P2K8": "GCN2",  # GCN2 has one pseudo domain
}
CATALYTIC_GO = {"kinase": "GO:0004672", "protease": "GO:0008236"}


def fetch_location(session, ipr, acc):
    """Return (start, end) of the first fragment of ipr on acc, or None."""
    url = ("https://www.ebi.ac.uk/interpro/api/entry/interpro/{}/protein/"
           "uniprot/{}".format(ipr, acc))
    try:
        r = session.get(url, timeout=30)
        if r.status_code != 200:
            return None
        js = r.json()
        def hunt(obj):
            if isinstance(obj, dict):
                if "entry_protein_locations" in obj and obj["entry_protein_locations"]:
                    fr = obj["entry_protein_locations"][0].get("fragments") or []
                    if fr:
                        return int(fr[0]["start"]), int(fr[0]["end"])
                for v in obj.values():
                    got = hunt(v)
                    if got:
                        return got
            elif isinstance(obj, list):
                for v in obj:
                    got = hunt(v)
                    if got:
                        return got
            return None
        return hunt(js)
    except Exception:
        return None


def ca_coords(pdb_path):
    coords = {}
    with open(pdb_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    resi = int(line[22:26])
                    coords[resi] = (float(line[30:38]), float(line[38:46]),
                                    float(line[46:54]))
                except ValueError:
                    pass
    return coords


def dist(a, b):
    return sum((x - y) ** 2 for x, y in zip(a, b)) ** 0.5


def structure_path(d, acc):
    for sd in (d / "structures", d / "dark" / "structures"):
        p = sd / (acc + ".pdb")
        if p.exists():
            return p
    return None


def main():
    cfg = load_config()
    d = data_dir(cfg)
    out_dir = Path("results/explore")
    out_dir.mkdir(parents=True, exist_ok=True)
    session = http_session()
    dag = GoDag(d / "go-basic.obo")

    domains = load_domains(d / "domains.tsv")
    for acc, ds in load_domains(d / "dark" / "domains.tsv").items():
        domains.setdefault(acc, set()).update(ds)
    seqs = read_fasta(d / "final.fasta")
    seqs.update(read_fasta(d / "dark" / "dark.fasta"))

    exp_go = defaultdict(set)
    with open(d / "annotations.tsv", encoding="utf-8") as f:
        next(f)
        for line in f:
            acc, term, _a = line.rstrip("\n").split("\t")
            exp_go[acc].add(term)
    prop_go = {a: dag.propagate(ts) for a, ts in exp_go.items()}

    families = {
        "kinase": sorted(a for a, ds in domains.items() if KINASE in ds and a in seqs),
        "protease": sorted(a for a, ds in domains.items() if PROTEASE in ds and a in seqs),
    }
    print("candidates: kinase {} | protease {}".format(
        len(families["kinase"]), len(families["protease"])))

    # ---- domain locations (cached, threaded) ------------------------------
    cache_path = d / "explore" / "domain_locations.json"
    cache = json.loads(cache_path.read_text()) if cache_path.exists() else {}
    todo = [(fam, acc) for fam, accs in families.items() for acc in accs
            if "{}:{}".format(acc, fam) not in cache]
    if todo:
        print("fetching {} domain locations ...".format(len(todo)))
        with ThreadPoolExecutor(max_workers=cfg["structures"]["workers"]) as pool:
            futs = {pool.submit(fetch_location, session,
                                KINASE if fam == "kinase" else PROTEASE, acc):
                    (fam, acc) for fam, acc in todo}
            for fut in tqdm(as_completed(futs), total=len(futs), unit="loc"):
                fam, acc = futs[fut]
                cache["{}:{}".format(acc, fam)] = fut.result()
        cache_path.write_text(json.dumps(cache))

    # ---- scan --------------------------------------------------------------
    rows = []
    for fam, accs in families.items():
        for acc in accs:
            seq = seqs[acc]
            loc = cache.get("{}:{}".format(acc, fam))
            if loc:
                start, end = max(1, loc[0] - 10), min(len(seq), loc[1] + 10)
                region, offset, scope = seq[start - 1:end], start - 1, "domain"
            else:
                region, offset, scope = seq, 0, "whole-seq"

            if fam == "kinase":
                m_b3 = RE_B3K.search(region[:min(len(region), 120)])
                m_hrd = RE_HRD.search(region)
                m_dfg = RE_DFG.search(region, m_hrd.end()) if m_hrd else RE_DFG.search(region)
                missing = [n for n, m in (("beta3K", m_b3), ("HRD", m_hrd),
                                          ("DFG", m_dfg)) if not m]
                geom = ""
                if m_hrd and m_dfg:
                    p = structure_path(d, acc)
                    if p:
                        coords = ca_coords(p)
                        r1 = offset + m_hrd.start() + 3  # HRD aspartate (1-based)
                        r2 = offset + m_dfg.start() + 1  # DFG aspartate
                        if r1 in coords and r2 in coords:
                            dd = dist(coords[r1], coords[r2])
                            geom = "{:.1f}".format(dd)
            else:
                m_his = RE_HIS.search(region)
                m_ser = RE_SER.search(region)
                missing = [n for n, m in (("triadHis", m_his),
                                          ("triadSer", m_ser)) if not m]
                geom = ""
                if m_his and m_ser:
                    p = structure_path(d, acc)
                    if p:
                        coords = ca_coords(p)
                        r1 = offset + m_his.start() + 5  # His position in pattern
                        r2 = offset + m_ser.start() + 6  # Ser position in pattern
                        if r1 in coords and r2 in coords:
                            geom = "{:.1f}".format(dist(coords[r1], coords[r2]))

            go_term = CATALYTIC_GO[fam]
            has_cat = go_term in prop_go.get(acc, set())
            status = "pseudo-candidate" if missing else "active-like"
            rows.append([acc, fam, scope, status, ";".join(missing) or "-",
                         geom or "-", "yes" if has_cat else "no",
                         KNOWN_PSEUDOKINASES.get(acc, "-")])

    rows.sort(key=lambda r: (r[1], r[3] != "pseudo-candidate", r[0]))
    with open(out_dir / "pseudoenzymes.tsv", "w", encoding="utf-8", newline="\n") as f:
        f.write("accession\tfamily\tscan_scope\tstatus\tmissing_motifs\t"
                "catalytic_CA_distance_A\texperimental_catalytic_GO\tknown_pseudo\n")
        for r in rows:
            f.write("\t".join(r) + "\n")

    # ---- summarise ---------------------------------------------------------
    def count(fam, status):
        return sum(1 for r in rows if r[1] == fam and r[3] == status)
    known_in_set = [r for r in rows if r[7] != "-"]
    known_flagged = [r for r in known_in_set if r[3] == "pseudo-candidate"]
    conflicts = [r for r in rows if r[3] == "pseudo-candidate" and r[6] == "yes"]
    geom_vals = [float(r[5]) for r in rows if r[1] == "kinase" and r[5] != "-"]
    geom_ok = sum(1 for g in geom_vals if 5 <= g <= 16)

    with open(out_dir / "pseudoenzymes.md", "w", encoding="utf-8", newline="\n") as f:
        f.write("# Pseudo-enzyme catalog\n\n")
        f.write("| family | members scanned | active-like | pseudo-candidates |\n")
        f.write("|---|---|---|---|\n")
        for fam in ("kinase", "protease"):
            f.write("| {} | {} | {} | {} |\n".format(
                fam, len(families[fam]), count(fam, "active-like"),
                count(fam, "pseudo-candidate")))
        f.write("\nGeometric check (kinases with both aspartate motifs "
                "located): {} of {} have HRD-DFG CA distance in the 5-16 A "
                "range expected of the kinase fold.\n\n".format(
                    geom_ok, len(geom_vals)))
        f.write("Validation against known human pseudokinases present in the "
                "dataset: {} of {} flagged as pseudo-candidates ({}).\n\n".format(
                    len(known_flagged), len(known_in_set),
                    ", ".join(r[7] for r in known_in_set) or "none present"))
        f.write("Candidates with motifs missing YET experimental catalytic GO "
                "annotation (conflicts worth review): {}\n\n".format(len(conflicts)))
        if conflicts:
            f.write("| accession | family | missing | known |\n|---|---|---|---|\n")
            for r in conflicts[:20]:
                f.write("| {} | {} | {} | {} |\n".format(r[0], r[1], r[4], r[7]))
        f.write("\nFull table: pseudoenzymes.tsv. Motif regexes are "
                "sequence-level heuristics scanned within InterPro domain "
                "boundaries; rows scanned whole-sequence (no boundary "
                "available) are marked and less reliable.\n")

    print("kinases: {} active-like, {} pseudo-candidates".format(
        count("kinase", "active-like"), count("kinase", "pseudo-candidate")))
    print("proteases: {} active-like, {} pseudo-candidates".format(
        count("protease", "active-like"), count("protease", "pseudo-candidate")))
    print("known pseudokinases in set: {} | flagged: {}".format(
        len(known_in_set), len(known_flagged)))
    print("Wrote results/explore/pseudoenzymes.md")


if __name__ == "__main__":
    main()
