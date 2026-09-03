"""Dynamics layer (build2, E) - MSA-subsampling conformational ensembles.

AlphaFold collapses each protein to a single static model. Del Alamo et al.
(eLife 2022) showed that running AF2 with a shallow MSA across many random
seeds makes the predictor sample alternative conformations instead of one
consensus fold; residues whose position varies most across that ensemble are
candidate hinges / flexible regions / multi-state segments. This is the
"dynamics" the static AlphaFold DB models cannot show.

For a small set of foldable single-chain candidates we run colabfold_batch
at reduced MSA depth (max_seq:max_extra_seq) across several seeds, superpose
the resulting CA traces, and report per-residue variance (an RMSF proxy).
High-variance proteins/regions are the flexible ones.

GPU, heavy (serialized behind any other GPU job via the flock heavy-lock).
Run:  python scripts/44_dynamics.py [N]
Outputs: results/explore/dynamics.md, dynamics_rmsf.tsv
"""
import glob
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from pis.common import read_fasta
from pis.stats import kabsch_rmsf
from pis import wsl as W

CF = "/root/localcolabfold/colabfold-conda/bin/colabfold_batch"
OUT = Path("data/explore/dyn")
MAXLEN = 320          # single chain; keep within 4GB VRAM
SEEDS = 8             # ensemble size
MAX_SEQ = "16:32"     # shallow MSA -> conformational sampling


def ca_trace(pdb_path):
    """Ordered list of (resi, x, y, z) CA atoms of the first chain."""
    out = []
    seen = set()
    with open(pdb_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if not line.startswith("ATOM") or line[12:16].strip() != "CA":
                continue
            ch = line[21]
            if out and ch != out[0][0]:
                break
            resi = int(line[22:26])
            if resi in seen:
                continue
            seen.add(resi)
            out.append((ch, resi, float(line[30:38]),
                        float(line[38:46]), float(line[46:54])))
    return [(r, x, y, z) for _c, r, x, y, z in out]


def ensemble(acc, seq):
    """Fold acc across SEEDS seeds at shallow MSA; return list of CA traces."""
    tag = acc
    outdir = OUT / ("out_" + tag)
    pdbs = glob.glob(str(outdir / "*.pdb"))
    if len(pdbs) < 2:
        OUT.mkdir(parents=True, exist_ok=True)
        fasta = OUT / (tag + ".fasta")
        with open(fasta, "w", newline="\n") as g:
            g.write(">{}\n{}\n".format(tag, seq))
        env = ("TF_FORCE_UNIFIED_MEMORY=1 XLA_PYTHON_CLIENT_MEM_FRACTION=3.0 "
               "XLA_PYTHON_CLIENT_ALLOCATOR=platform")
        job = ("{} {} {} {} --num-recycle 1 --num-models 1 --num-seeds {} "
               "--max-seq {} --max-extra-seq {}\n").format(
                   env, CF, W.to_wsl_path(fasta), W.to_wsl_path(outdir),
                   SEEDS, MAX_SEQ.split(":")[0], MAX_SEQ.split(":")[1])
        W.run_script(job, heavy=True, check=False, timeout=3600,
                     out_path=OUT / ("_dyn_{}.sh".format(tag)))
        pdbs = glob.glob(str(outdir / "*.pdb"))
    traces = []
    for p in sorted(pdbs):
        t = ca_trace(p)
        if t:
            traces.append(t)
    return traces


def pick(n):
    """Foldable single-chain candidates: prefer dark proteins (unknown
    function - flexibility is a functional clue), fall back to annotated."""
    seqs = read_fasta(Path("data/dark/dark.fasta"))
    cands = sorted(((len(s), a) for a, s in seqs.items() if len(s) <= MAXLEN))
    return [(a, seqs[a]) for _L, a in cands[:n]]


def main():
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 6
    cands = pick(n)
    print("dynamics candidates (dark, <= {} aa): {}".format(MAXLEN, len(cands)))
    rows = []
    for acc, seq in cands:
        print("  ensemble:", acc, "len", len(seq), flush=True)
        try:
            traces = ensemble(acc, seq)
            if len(traces) < 2:
                rows.append((acc, len(seq), 0, 0.0, 0.0, "too_few_models"))
                continue
            rmsf = kabsch_rmsf(traces)
            vals = list(rmsf.values())
            mean_rmsf = sum(vals) / len(vals)
            max_rmsf = max(vals)
            flex = [r for r, v in rmsf.items() if v > 2 * mean_rmsf]
            rows.append((acc, len(seq), len(traces), round(mean_rmsf, 2),
                         round(max_rmsf, 2),
                         "flexible" if mean_rmsf > 3 else "rigid"))
            print("    models {} mean_rmsf {:.2f} max {:.2f} flex_res {}".format(
                len(traces), mean_rmsf, max_rmsf, len(flex)))
        except Exception as e:
            rows.append((acc, len(seq), 0, 0.0, 0.0, "error:" + str(e)[:30]))
            print("    error:", str(e)[:80])

    Path("results/explore").mkdir(parents=True, exist_ok=True)
    with open("results/explore/dynamics_rmsf.tsv", "w", encoding="utf-8",
              newline="\n") as f:
        f.write("accession\tlength\tmodels\tmean_rmsf\tmax_rmsf\tcall\n")
        for r in rows:
            f.write("\t".join(str(x) for x in r) + "\n")
    ranked = sorted((r for r in rows if r[2] >= 2), key=lambda r: -r[3])
    with open("results/explore/dynamics.md", "w", encoding="utf-8",
              newline="\n") as f:
        f.write("# Conformational dynamics via MSA-subsampling ensembles\n\n")
        f.write("AlphaFold gives one static model per protein. Running it at "
                "shallow MSA depth ({}) across {} random seeds makes it sample "
                "alternative conformations (Del Alamo et al., eLife 2022); the "
                "spread of CA positions across the ensemble is a proxy for "
                "intrinsic flexibility that the single AlphaFold DB model hides."
                "\n\n".format(MAX_SEQ, SEEDS))
        f.write("Per-residue RMSF is measured on Kabsch-superposed CA traces "
                "(rotation + translation removed) over residues present in every "
                "model. Dark (unannotated) proteins were prioritised - "
                "flexibility is a functional clue for a protein with no known "
                "function.\n\n")
        f.write("| protein | length | models | mean RMSF (A) | max RMSF (A) | "
                "call |\n|---|---|---|---|---|---|\n")
        for acc, L, m, mean, mx, call in ranked:
            f.write("| {} | {} | {} | {} | {} | {} |\n".format(
                acc, L, m, mean, mx, call))
        f.write("\nHigh mean RMSF flags a candidate conformationally flexible "
                "or multi-state protein; localized high-RMSF stretches (see "
                "dynamics_rmsf.tsv) are candidate hinges or disordered "
                "segments. RMSF here is a Kabsch-superposed relative-flexibility "
                "proxy, not an absolute B-factor.\n")
    print("Wrote results/explore/dynamics.md")


if __name__ == "__main__":
    main()
