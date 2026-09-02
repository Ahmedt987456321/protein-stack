"""CAFA-style protein-centric Fmax.

Fast path: each protein's predictions are bucketed once into per-threshold
cumulative counts (predicted / correct), so sweeping the 100 thresholds costs
O(1) per protein per threshold instead of rebuilding term sets.
"""

N_T = 100
THRESHOLDS = [i / N_T for i in range(1, N_T + 1)]


def protein_curve(pred, truth):
    """pred: {term: score}, truth: set of terms.

    Returns (npred, ncorr): for threshold index j (= threshold (j+1)/100),
    npred[j] = #terms predicted with score >= threshold, ncorr[j] of them true.
    """
    npred = [0] * N_T
    ncorr = [0] * N_T
    for term, s in pred.items():
        k = int(s * N_T + 1e-9)  # number of thresholds satisfied by score s
        if k <= 0:
            continue
        if k > N_T:
            k = N_T
        npred[k - 1] += 1
        if term in truth:
            ncorr[k - 1] += 1
    for j in range(N_T - 2, -1, -1):
        npred[j] += npred[j + 1]
        ncorr[j] += ncorr[j + 1]
    return npred, ncorr


def fmax_from_curves(curves, truth_sizes):
    """curves: list of (npred, ncorr); truth_sizes: per-protein |truth|.

    Precision at t is averaged over proteins with >=1 prediction at t;
    recall at t is averaged over all proteins. Returns (fmax, best_threshold).
    """
    n = len(curves)
    if n == 0:
        return 0.0, 0.0
    best_f, best_t = 0.0, 0.0
    for j in range(N_T):
        prec_sum, prec_n, rec_sum = 0.0, 0, 0.0
        for (npred, ncorr), ts in zip(curves, truth_sizes):
            if npred[j] > 0:
                prec_sum += ncorr[j] / npred[j]
                prec_n += 1
            rec_sum += ncorr[j] / ts
        if prec_n == 0:
            continue
        pr = prec_sum / prec_n
        rc = rec_sum / n
        if pr + rc > 0:
            f = 2 * pr * rc / (pr + rc)
            if f > best_f:
                best_f, best_t = f, THRESHOLDS[j]
    return best_f, best_t


def fmax(proteins, preds, truth):
    """Drop-in equivalent of the v0.1 implementation."""
    curves, sizes = [], []
    for acc in proteins:
        curves.append(protein_curve(preds.get(acc, {}), truth[acc]))
        sizes.append(len(truth[acc]))
    return fmax_from_curves(curves, sizes)
