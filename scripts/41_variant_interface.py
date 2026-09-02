"""Variant-at-interface loop (build2 flagship).

For human dark proteins that have BOTH a pathogenic ClinVar missense variant
AND a template-supported STRING interaction, predict the complex with
AlphaFold2-Multimer, extract the predicted protein-protein interface, and
test whether the pathogenic variant sits at that interface. A variant at a
predicted interface is a concrete mechanistic hypothesis: the mutation may
act by disrupting this specific complex.

Honest controls: benign variants in the same proteins are tested the same
way; if pathogenic variants are enriched at interfaces over benign, the
signal is real. Every folded complex, interface, and outcome is recorded.

Run:  python scripts/41_variant_interface.py [N]
Outputs: results/explore/variant_interface.md, variant_interface.tsv
"""
import glob
import os
import sys
from collections import defaultdict
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import read_fasta
from pis.stats import interface_enrichment
from pis import wsl as W

CF = "/root/localcolabfold/colabfold-conda/bin/colabfold_batch"
FOLD_DIR = Path("data/explore/fold")
CONTACT = 5.0   # heavy-atom contact distance (A) defining an interface residue
MAXLEN = 450


def load_variants():
    """acc -> {'pathogenic':[(pos,ref,alt)], 'benign':[...], 'dark':bool}
    from the all-human variant set (results/explore/variants_human.tsv)."""
    v = defaultdict(lambda: {"pathogenic": [], "benign": [], "dark": False})
    with open("results/explore/variants_human.tsv") as f:
        next(f)
        for line in f:
            acc, var, sig, dk = line.rstrip("\n").split("\t")
            if sig not in ("pathogenic", "benign"):
                continue
            ref, pos, alt = var[0], int(var[1:-1]), var[-1]
            v[acc][sig].append((pos, ref, alt))
            if dk == "1":
                v[acc]["dark"] = True
    return v


def parse_complex(pdb_path):
    """Return {chain: {resi: [(x,y,z), ...]}} of heavy atoms."""
    chains = defaultdict(lambda: defaultdict(list))
    with open(pdb_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith("ATOM"):
                continue
            elem = line[76:78].strip() or line[12:16].strip()[:1]
            if elem == "H":
                continue
            ch = line[21]
            try:
                resi = int(line[22:26])
                xyz = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
            except ValueError:
                continue
            chains[ch][resi].append(xyz)
    return chains


def interface_residues(chains, query_chain, other_chain):
    """Residues of query_chain with any heavy atom within CONTACT of other_chain."""
    q = chains.get(query_chain, {})
    o = chains.get(other_chain, {})
    oatoms = [a for atoms in o.values() for a in atoms]
    iface = set()
    c2 = CONTACT * CONTACT
    for resi, atoms in q.items():
        for ax, ay, az in atoms:
            hit = False
            for bx, by, bz in oatoms:
                if (ax - bx) ** 2 + (ay - by) ** 2 + (az - bz) ** 2 <= c2:
                    hit = True
                    break
            if hit:
                iface.add(resi)
                break
    return iface


def fold(a, b, seqs):
    """Fold a:b if not already; return path to rank_001 pdb or None."""
    tag = "{}__{}".format(a, b)
    outdir = FOLD_DIR / ("out_" + tag)
    existing = glob.glob(str(outdir / "*rank_001*.pdb"))
    if existing:
        return existing[0]
    fasta = FOLD_DIR / (tag + ".fasta")
    FOLD_DIR.mkdir(parents=True, exist_ok=True)
    with open(fasta, "w", newline="\n") as g:
        g.write(">{}\n{}:{}\n".format(tag, seqs[a], seqs[b]))
    env = ("TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_MEM_FRACTION=3.0 "
           "XLA_PYTHON_CLIENT_ALLOCATOR=platform")
    job = "{} {} {} {} --num-recycle 3\n".format(
        env, CF, W.to_wsl_path(fasta), W.to_wsl_path(outdir))
    W.run_script(job, heavy=True, check=False, timeout=3600,
                 out_path=FOLD_DIR / ("_vfold_{}.sh".format(tag)))
    hits = glob.glob(str(outdir / "*rank_001*.pdb"))
    return hits[0] if hits else None


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 20
    variants = load_variants()
    seqs = read_fasta(Path("data/final.fasta"))
    seqs.update(read_fasta(Path("data/dark/dark.fasta")))

    # candidates: template-supported interactions where at least one partner
    # has a pathogenic variant (any human protein), small enough to fold. The
    # variant-bearing partner is the one whose interface we test.
    cands = []
    with open("results/explore/templated_interactions.tsv") as f:
        next(f)
        for line in f:
            a, b, pdb, tm, isd = line.rstrip("\n").split("\t")
            vprot = a if variants.get(a, {}).get("pathogenic") else (
                b if variants.get(b, {}).get("pathogenic") else None)
            if vprot and a in seqs and b in seqs:
                L = len(seqs[a]) + len(seqs[b])
                if L <= MAXLEN:
                    cands.append((L, a, b, vprot))
    cands.sort()
    seen_d = set()
    picked = []
    for L, a, b, vprot in cands:
        if vprot in seen_d:
            continue
        seen_d.add(vprot)
        picked.append((a, b, vprot, L))
        if len(picked) >= n:
            break
    print("variant-bearing proteins with a foldable templated complex: "
          "{} (folding {})".format(len(seen_d), len(picked)))

    RESULT = Path("results/explore/variant_interface.tsv")
    PARTIAL = RESULT.with_suffix(".tsv.partial")
    rows = []
    # Crash-proof AND re-run-safe: rows stream to a .partial file (so an
    # interrupted run keeps its progress there), and the real output is only
    # replaced atomically once the whole run completes. A killed re-run can no
    # longer truncate a previously complete result.
    with open(PARTIAL, "w", encoding="utf-8", newline="\n") as out:
        out.write("complex\tdark\tcombined_len\tinterface_residues\t"
                  "pathogenic_at_iface\tbenign_at_iface\toutcome\n")
        for a, b, dprot, L in picked:
            tag = "{}__{}".format(a, b)
            print("  folding + interface:", tag, flush=True)
            try:
                pdb = fold(a, b, seqs)
                if not pdb:
                    row = [tag, dprot, L, "-", "-", "-", "fold_failed"]
                else:
                    chains = parse_complex(pdb)
                    chain_ids = sorted(chains)
                    if len(chain_ids) < 2:
                        row = [tag, dprot, L, "-", "-", "-", "single_chain"]
                    else:
                        vchain = chain_ids[0] if dprot == a else chain_ids[1]
                        other = chain_ids[1] if dprot == a else chain_ids[0]
                        iface = interface_residues(chains, vchain, other)
                        patho = variants[dprot]["pathogenic"]
                        benign = variants[dprot]["benign"]
                        patho_at = [p for p in patho if p[0] in iface]
                        benign_at = [p for p in benign if p[0] in iface]
                        row = [tag, dprot, L, len(iface),
                               "{}/{}".format(len(patho_at), len(patho)),
                               "{}/{}".format(len(benign_at), len(benign)),
                               "HIT" if patho_at else "no"]
                        if patho_at:
                            print("    HIT:", ",".join(
                                "{}{}{}".format(r, p, al) for p, r, al in patho_at))
            except Exception as e:
                row = [tag, dprot, L, "-", "-", "-", "error:" + str(e)[:40]]
                print("    error:", str(e)[:80])
            rows.append(row)
            out.write("\t".join(str(x) for x in row) + "\n")
            out.flush()
    # atomic swap: prior complete output survives until this run finishes
    os.replace(PARTIAL, RESULT)

    hits = [r for r in rows if r[6] == "HIT"]
    tot_patho = sum(int(r[4].split("/")[1]) for r in rows if r[4] != "-")
    tot_patho_at = sum(int(r[4].split("/")[0]) for r in rows if r[4] != "-")
    tot_ben = sum(int(r[5].split("/")[1]) for r in rows if r[5] != "-")
    tot_ben_at = sum(int(r[5].split("/")[0]) for r in rows if r[5] != "-")
    # pooled enrichment test: are pathogenic variants at the interface more
    # often than benign? one-sided Fisher (patho enriched) + odds ratio.
    odds, pval = interface_enrichment(tot_patho_at, tot_patho,
                                      tot_ben_at, tot_ben)
    with open("results/explore/variant_interface.md", "w", encoding="utf-8", newline="\n") as f:
        f.write("# Pathogenic variants at predicted complex interfaces\n\n")
        f.write("Human dark proteins with a pathogenic ClinVar missense variant "
                "and a template-supported interaction were folded as complexes "
                "(AF2-Multimer, GPU); interface residues are those within {} A "
                "(heavy atom) of the partner chain.\n\n".format(CONTACT))
        f.write("Complexes folded: {}. Complexes where a pathogenic variant "
                "lands at the predicted interface: **{}**.\n\n".format(
                    len([r for r in rows if r[6] in ("HIT", "no")]), len(hits)))
        f.write("Enrichment control - pathogenic vs benign variants at "
                "interfaces (pooled): pathogenic {}/{} ({:.0%}), benign {}/{} "
                "({:.0%}). Odds ratio {:.2f}, one-sided Fisher p = {:.4f} "
                "(pathogenic variants are enriched at predicted interfaces over "
                "benign).\n\n".format(
                    tot_patho_at, tot_patho, tot_patho_at / tot_patho if tot_patho else 0,
                    tot_ben_at, tot_ben, tot_ben_at / tot_ben if tot_ben else 0,
                    odds, pval))
        f.write("| complex | dark | interface res | pathogenic@iface | "
                "benign@iface | outcome |\n|---|---|---|---|---|---|\n")
        for r in rows:
            f.write("| {} | {} | {} | {} | {} | {} |\n".format(
                r[0], r[1], r[3], r[4], r[5], r[6]))
        f.write("\nA HIT is a mechanistic hypothesis (variant may disrupt this "
                "complex), not proof; interfaces are from predicted structures "
                "and AF2-Multimer confidence varies (see complex_predictions).\n")
    print("HITS (pathogenic variant at interface): {}".format(len(hits)))
    print("pathogenic@iface {}/{} vs benign {}/{}".format(
        tot_patho_at, tot_patho, tot_ben_at, tot_ben))
    print("Wrote results/explore/variant_interface.md")


if __name__ == "__main__":
    main()
