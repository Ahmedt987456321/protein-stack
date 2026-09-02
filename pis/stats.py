"""Load-bearing statistics for the published results.

These two functions produce headline numbers in the paper (the pathogenic-
variant interface enrichment and the conformational-ensemble RMSF), so they
live here in the tested library rather than inline in analysis scripts. Both
are pure and deterministic; see tests/test_core.py for hand-computed checks.
"""
from __future__ import annotations


def interface_enrichment(patho_hit, patho_tot, benign_hit, benign_tot):
    """One-sided Fisher exact test that pathogenic variants sit at a predicted
    interface more often than benign variants.

    Args are counts: (variants of each class at an interface, total of each
    class). Returns (odds_ratio, p_value) where p is the one-sided probability
    (alternative: pathogenic enriched). Returns (nan, nan) if either class has
    no variants, since the ratio is undefined.

    2x2 contingency table:
                        at interface     not at interface
        pathogenic      patho_hit        patho_tot - patho_hit
        benign          benign_hit       benign_tot - benign_hit
    """
    if patho_tot <= 0 or benign_tot <= 0:
        return float("nan"), float("nan")
    for n, name in ((patho_hit, "patho_hit"), (patho_tot, "patho_tot"),
                    (benign_hit, "benign_hit"), (benign_tot, "benign_tot")):
        if n < 0:
            raise ValueError("{} must be non-negative, got {}".format(name, n))
    if patho_hit > patho_tot or benign_hit > benign_tot:
        raise ValueError("hits cannot exceed totals")
    from scipy.stats import fisher_exact
    table = [[patho_hit, patho_tot - patho_hit],
             [benign_hit, benign_tot - benign_hit]]
    odds, p = fisher_exact(table, alternative="greater")
    return float(odds), float(p)


def kabsch_rmsf(traces):
    """Per-residue RMSF over an ensemble of CA traces after rigid-body
    superposition (Kabsch/SVD rotation + centroid translation) onto the first
    trace as reference.

    Rotational alignment is essential: models that fold identically but are
    written in different absolute orientations (as AlphaFold/ColabFold emits
    them across seeds) otherwise show large spurious flexibility. Only residues
    present in every trace are scored.

    Args:
        traces: list of traces; each trace is a list of (resi, x, y, z) tuples
            for the CA atoms of one model.

    Returns:
        dict {resi: rmsf_angstrom}. Empty dict if fewer than two traces or
        fewer than three residues are shared across all traces (SVD needs a
        non-degenerate point set).
    """
    import numpy as np
    if len(traces) < 2:
        return {}
    common = set(r for r, *_ in traces[0])
    for t in traces[1:]:
        common &= set(r for r, *_ in t)
    common = sorted(common)
    if len(common) < 3:
        return {}
    mats = []
    for t in traces:
        d = {r: (x, y, z) for r, x, y, z in t if r in common}
        mats.append(np.array([d[r] for r in common], dtype=float))
    ref = mats[0] - mats[0].mean(axis=0)
    aligned = []
    for m in mats:
        mc = m - m.mean(axis=0)
        # Kabsch rotation minimizing RMSD of mc onto ref, with reflection guard
        H = mc.T @ ref
        U, _S, Vt = np.linalg.svd(H)
        d = np.sign(np.linalg.det(Vt.T @ U.T))
        D = np.diag([1.0, 1.0, d])
        R = Vt.T @ D @ U.T
        aligned.append(mc @ R.T)
    stack = np.stack(aligned)                     # models x N x 3
    mean_xyz = stack.mean(axis=0)                 # N x 3
    dev = ((stack - mean_xyz) ** 2).sum(axis=2).mean(axis=0)  # N
    rmsf = np.sqrt(dev)
    return {r: float(v) for r, v in zip(common, rmsf)}
