"""Batch AlphaFold2-Multimer prediction of dark-protein interactions on GPU.

Selects the smallest, most diverse dark-involving STRING pairs (one pair per
dark protein, combined length <= MAXLEN so they fit a 4GB card), folds each
with ColabFold on the GPU, and records ipTM/pTM. Bulky outputs are pruned to
the top-ranked model to save disk. Failures (OOM, errors) are recorded, not
dropped.

Run:  python scripts/40_batch_fold.py [N]
"""
import glob
import json
import os
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import read_fasta
from pis import wsl as W

CF = "/root/localcolabfold/colabfold-conda/bin/colabfold_batch"
MAXLEN = 450
FOLD_DIR = Path("data/explore/fold")
RESULT = Path("results/explore/complex_predictions.tsv")


def pick_candidates(n):
    """Prefer dark interactions with a PDB structural template (strongest
    candidates: STRING + template agree), smallest first, one per dark
    protein. Falls back to plain smallest dark STRING pairs if the template
    file is absent."""
    seqs = read_fasta(Path("data/final.fasta"))
    seqs.update(read_fasta(Path("data/dark/dark.fasta")))
    dark = set(Path("data/dark/accessions.txt").read_text().split())
    have = set(seqs)
    templ = Path("results/explore/templated_interactions.tsv")
    rows = []
    if templ.exists():
        with open(templ) as f:
            next(f)
            for line in f:
                a, b, pdb, tm, isd = line.rstrip("\n").split("\t")
                if isd == "1" and a in have and b in have:
                    L = len(seqs[a]) + len(seqs[b])
                    if L <= MAXLEN:
                        rows.append((L, a, b, pdb))
    if not rows:  # fallback: plain STRING dark pairs
        with open("data/string_edges.tsv") as f:
            next(f)
            for line in f:
                a, b, s = line.rstrip("\n").split("\t")
                if a in have and b in have and (a in dark or b in dark):
                    L = len(seqs[a]) + len(seqs[b])
                    if L <= MAXLEN:
                        rows.append((L, a, b, "-"))
    rows.sort()
    picked, used_dark = [], set()
    for L, a, b, pdb in rows:
        dprot = a if a in dark else b
        if dprot in used_dark:
            continue
        used_dark.add(dprot)
        picked.append((a, b, L, pdb, dprot))
        if len(picked) >= n:
            break
    return picked, seqs


def read_scores(outdir):
    for pat in ("*scores_rank_001*.json", "*rank_001*.json"):
        hits = sorted(glob.glob(os.path.join(outdir, pat)))
        if hits:
            j = json.load(open(hits[0]))
            plddt = j.get("plddt", [])
            return {"iptm": j.get("iptm"), "ptm": j.get("ptm"),
                    "plddt_mean": round(sum(plddt) / len(plddt), 1) if plddt else None}
    return None


def prune(outdir):
    """Keep only rank_001 model + its scores; delete the rest to save disk."""
    for f in glob.glob(os.path.join(outdir, "*")):
        if "rank_001" in os.path.basename(f) or f.endswith(".fasta"):
            continue
        try:
            os.remove(f)
        except OSError:
            pass


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 15
    FOLD_DIR.mkdir(parents=True, exist_ok=True)
    RESULT.parent.mkdir(parents=True, exist_ok=True)
    picked, seqs = pick_candidates(n)
    print("folding {} candidates (combined length <= {}):".format(len(picked), MAXLEN))

    done = {}
    if RESULT.exists():
        with open(RESULT) as f:
            next(f, None)
            for line in f:
                done[line.split("\t")[0]] = True

    with open(RESULT, "a", encoding="utf-8", newline="\n") as out:
        if out.tell() == 0:
            out.write("pair\tprotein_a\tprotein_b\tdark\tcombined_len\t"
                      "template_pdb\tiptm\tptm\tplddt_mean\tstatus\n")
        for a, b, L, sc, dprot in picked:
            tag = "{}__{}".format(a, b)
            if tag in done:
                print("  skip (done):", tag)
                continue
            fasta = FOLD_DIR / (tag + ".fasta")
            with open(fasta, "w", newline="\n") as g:
                g.write(">{}\n{}:{}\n".format(tag, seqs[a], seqs[b]))
            outdir = FOLD_DIR / ("out_" + tag)
            env = ("TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_MEM_FRACTION=3.0 "
                   "XLA_PYTHON_CLIENT_ALLOCATOR=platform")
            job = "{} {} {} {} --num-recycle 3\n".format(
                env, CF, W.to_wsl_path(fasta), W.to_wsl_path(outdir))
            print("  folding", tag, "(len {})".format(L), "...", flush=True)
            # heavy=True serialises with any concurrent Foldseek search (OOM
            # guard); routed through pis.wsl so paths are never mangled.
            r = W.run_script(job, heavy=True, check=False, timeout=3600,
                             out_path=FOLD_DIR / ("_fold_{}.sh".format(tag)))
            sc_d = read_scores(str(outdir)) if outdir.exists() else None
            if sc_d and sc_d.get("iptm") is not None:
                status = "ok"
                out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(
                    tag, a, b, dprot, L, sc, sc_d["iptm"], sc_d["ptm"],
                    sc_d["plddt_mean"], status))
                print("    ipTM={} pTM={}".format(sc_d["iptm"], sc_d["ptm"]))
                prune(str(outdir))
            else:
                tail = (r.stderr or r.stdout or "")[-300:].replace("\n", " ")
                status = "OOM" if "out of memory" in tail.lower() or "resource_exhausted" in tail.lower() else "fail"
                out.write("{}\t{}\t{}\t{}\t{}\t{}\t-\t-\t-\t{}\n".format(
                    tag, a, b, dprot, L, sc, status))
                print("    {}: {}".format(status, tail[-160:]))
            out.flush()
    print("Wrote", RESULT)


if __name__ == "__main__":
    main()
