"""Exploration 5 - predicted disorder versus annotation status.

Compares AlphaFold confidence between experimentally annotated and dark
proteins, without survivor bias: the dark models that failed the pLDDT gate
in step 11 are re-derived (same deterministic sample), fetched transiently,
measured, and deleted.

Metrics: mean pLDDT (pre-gate, both sets) and fraction of residues with
pLDDT < 50 (a standard disorder proxy; post-gate kept sets, both on disk).
Significance: rank-based test via 10,000-fold permutation of group labels
on the median difference.

Outputs: results/explore/disorder.md, results/explore/disorder_stats.csv
"""
import csv
import random
import statistics
import sys
import tempfile
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from tqdm import tqdm

from pis.common import data_dir, http_session, load_config
from pis.fetch import UNIPROT_BATCH, fetch_alphafold, fetch_uniprot_batch
from pis.go import parse_gaf_full


def plddt_stats(pdb_path):
    vals = []
    with open(pdb_path, encoding="utf-8", errors="replace") as f:
        for line in f:
            if line.startswith("ATOM") and line[12:16].strip() == "CA":
                try:
                    vals.append(float(line[60:66]))
                except ValueError:
                    pass
    if not vals:
        return None
    return (sum(vals) / len(vals), sum(v < 50 for v in vals) / len(vals))


def perm_test(a, b, rng, n=10000):
    obs = statistics.median(a) - statistics.median(b)
    pooled = a + b
    na = len(a)
    ge = 0
    for _ in range(n):
        rng.shuffle(pooled)
        if abs(statistics.median(pooled[:na]) - statistics.median(pooled[na:])) >= abs(obs):
            ge += 1
    return obs, (ge + 1) / (n + 1)


def main():
    cfg = load_config()
    d = data_dir(cfg)
    out_dir = Path("results/explore")
    out_dir.mkdir(parents=True, exist_ok=True)
    session = http_session()
    rng = random.Random(cfg["seed"] + 30)

    # ---- reconstruct the step-11 dark sample deterministically ------------
    exp_codes = set(cfg["evidence_codes"])
    has_exp = set()
    nonexp = set()
    for sp in cfg["species"]:
        gaf = d / "gaf" / "goa_{}.gaf.gz".format(sp["name"])
        for acc, _t, _a, ev in parse_gaf_full(gaf):
            (has_exp if ev in exp_codes else nonexp).add(acc)
    dark_all = sorted(nonexp - has_exp)
    rng11 = random.Random(cfg["seed"] + 3)
    rng11.shuffle(dark_all)
    sample = dark_all[: cfg["dark"]["oversample"]]

    print("Refetching sequences to reproduce the length filter ...")
    seqs = {}
    for i in tqdm(range(0, len(sample), UNIPROT_BATCH), unit="batch"):
        seqs.update(fetch_uniprot_batch(session, sample[i:i + UNIPROT_BATCH]))
    lo, hi = cfg["sampling"]["min_len"], cfg["sampling"]["max_len"]
    length_ok = sorted(a for a, s in seqs.items() if lo <= len(s) <= hi)
    kept = set((d / "dark" / "accessions.txt").read_text(encoding="utf-8").split())
    failed = [a for a in length_ok if a not in kept]
    print("length-passed {} | kept {} | gate-failed to refetch {}".format(
        len(length_ok), len(kept), len(failed)))

    # ---- transient fetch of gate-failed dark models -----------------------
    dark_fail_means = []
    with tempfile.TemporaryDirectory() as td:
        tdir = Path(td)
        def one(acc):
            acc2, plddt, ok = fetch_alphafold(session, acc, tdir, min_plddt=-1.0)
            f = tdir / (acc + ".pdb")
            st = plddt_stats(f) if f.exists() else None
            f.unlink(missing_ok=True)
            return st[0] if st else None
        with ThreadPoolExecutor(max_workers=cfg["structures"]["workers"]) as pool:
            futs = [pool.submit(one, a) for a in failed]
            for fut in tqdm(as_completed(futs), total=len(futs), unit="model"):
                m = fut.result()
                if m is not None:
                    dark_fail_means.append(m)
    print("gate-failed models measured:", len(dark_fail_means))

    # ---- assemble distributions ------------------------------------------
    ann_means_all = []
    with open(d / "structures.csv", encoding="utf-8") as f:
        next(f)
        for row in csv.reader(f):
            m = float(row[1])
            if m >= 0:
                ann_means_all.append(m)
    dark_kept_stats = []
    for p in (d / "dark" / "structures").glob("*.pdb"):
        st = plddt_stats(p)
        if st:
            dark_kept_stats.append(st)
    ann_kept_stats = []
    for p in (d / "structures").glob("*.pdb"):
        st = plddt_stats(p)
        if st:
            ann_kept_stats.append(st)
    dark_means_all = [s[0] for s in dark_kept_stats] + dark_fail_means

    obs_pre, p_pre = perm_test(ann_means_all[:], dark_means_all[:], rng)
    ann_d50 = [s[1] for s in ann_kept_stats]
    dark_d50 = [s[1] for s in dark_kept_stats]
    obs_d50, p_d50 = perm_test(ann_d50[:], dark_d50[:], rng)

    def med(x): return statistics.median(x)
    rows = [
        ["pre-gate mean pLDDT", "annotated", len(ann_means_all), round(med(ann_means_all), 2)],
        ["pre-gate mean pLDDT", "dark", len(dark_means_all), round(med(dark_means_all), 2)],
        ["kept fraction pLDDT<50", "annotated", len(ann_d50), round(med(ann_d50), 4)],
        ["kept fraction pLDDT<50", "dark", len(dark_d50), round(med(dark_d50), 4)],
    ]
    with open(out_dir / "disorder_stats.csv", "w", encoding="utf-8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["metric", "group", "n", "median"])
        w.writerows(rows)
        w.writerow(["pre-gate median diff (ann-dark)", "", "", round(obs_pre, 2)])
        w.writerow(["pre-gate permutation p", "", "", p_pre])
        w.writerow(["kept d50 median diff (ann-dark)", "", "", round(obs_d50, 4)])
        w.writerow(["kept d50 permutation p", "", "", p_d50])

    with open(out_dir / "disorder.md", "w", encoding="utf-8", newline="\n") as f:
        f.write("# Predicted disorder vs annotation status\n\n")
        f.write("Pre-gate mean pLDDT (no survivor bias; dark gate failures "
                "re-measured transiently):\n\n")
        f.write("| group | n | median mean-pLDDT |\n|---|---|---|\n")
        f.write("| annotated | {} | {:.2f} |\n".format(len(ann_means_all), med(ann_means_all)))
        f.write("| dark | {} | {:.2f} |\n\n".format(len(dark_means_all), med(dark_means_all)))
        f.write("Median difference {:.2f}, permutation p = {}.\n\n".format(obs_pre, p_pre))
        f.write("Post-gate disorder proxy (fraction of residues with pLDDT "
                "< 50, kept models only):\n\n")
        f.write("| group | n | median fraction |\n|---|---|---|\n")
        f.write("| annotated | {} | {:.4f} |\n".format(len(ann_d50), med(ann_d50)))
        f.write("| dark | {} | {:.4f} |\n\n".format(len(dark_d50), med(dark_d50)))
        f.write("Median difference {:.4f}, permutation p = {}.\n".format(obs_d50, p_d50))
    print("pre-gate: ann {:.2f} vs dark {:.2f} (p={})".format(
        med(ann_means_all), med(dark_means_all), p_pre))
    print("kept d50: ann {:.4f} vs dark {:.4f} (p={})".format(
        med(ann_d50), med(dark_d50), p_d50))
    print("Wrote results/explore/disorder.md")


if __name__ == "__main__":
    main()
