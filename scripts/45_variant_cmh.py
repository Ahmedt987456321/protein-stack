"""Cluster-aware re-analysis of the variant-at-interface result (review point).

A plain Fisher test pools all variants and treats them as independent, but
variants from the same protein/complex are correlated. This stratifies by
complex: a Cochran-Mantel-Haenszel test (common odds ratio + CMH chi-square)
combines the per-complex 2x2 tables, and a within-complex permutation gives a
distribution-free check. Both respect the clustering the reviewer flagged.

Reads results/explore/variant_interface.tsv; writes results/explore/variant_cmh.md.
"""
import math
import random
from pathlib import Path


def load():
    """Per complex: (patho_hit, patho_tot, benign_hit, benign_tot)."""
    rows = []
    with open("results/explore/variant_interface.tsv") as f:
        next(f)
        for line in f:
            c = line.rstrip("\n").split("\t")
            if len(c) < 6 or c[4] == "-":
                continue
            ph, pt = (int(x) for x in c[4].split("/"))
            bh, bt = (int(x) for x in c[5].split("/"))
            rows.append((ph, pt, bh, bt))
    return rows


def mantel_haenszel(strata):
    """MH common odds ratio + CMH chi-square (1 df) with continuity correction.
    Each stratum: a=patho@iface, b=patho not, c=benign@iface, d=benign not."""
    num_or = den_or = 0.0
    sum_a = sum_ea = sum_va = 0.0
    used = 0
    for ph, pt, bh, bt in strata:
        a, b, c, d = ph, pt - ph, bh, bt - bh
        n = a + b + c + d
        if n == 0:
            continue
        r1, r2 = a + b, c + d          # pathogenic, benign totals
        c1 = a + c                     # at-interface total
        if r1 == 0 or r2 == 0 or c1 == 0 or c1 == n:
            continue                   # stratum carries no information
        used += 1
        num_or += a * d / n
        den_or += b * c / n
        sum_a += a
        sum_ea += r1 * c1 / n
        sum_va += r1 * r2 * c1 * (n - c1) / (n * n * (n - 1)) if n > 1 else 0.0
    mh_or = num_or / den_or if den_or else float("nan")
    chi = (abs(sum_a - sum_ea) - 0.5) ** 2 / sum_va if sum_va else 0.0
    # one-sided p (pathogenic enriched): normal approx on the signed z
    z = (sum_a - sum_ea) / math.sqrt(sum_va) if sum_va else 0.0
    p_two = math.erfc(math.sqrt(chi / 2)) if chi > 0 else 1.0
    p_one = p_two / 2 if z > 0 else 1 - p_two / 2
    return mh_or, chi, p_one, used


def perm_test(strata, n=20000, seed=0):
    """Within-complex permutation: in each stratum redistribute the interface
    hits among that stratum's variants (hypergeometric), keeping margins. Test
    statistic = pooled (patho@iface rate - benign@iface rate)."""
    rng = random.Random(seed)

    def stat(assign):
        pa = pt_ = ba = bt_ = 0
        for (ph, pt, bh, bt), (a, c) in zip(strata, assign):
            pa += a; pt_ += pt; ba += c; bt_ += bt
        rp = pa / pt_ if pt_ else 0.0
        rb = ba / bt_ if bt_ else 0.0
        return rp - rb

    obs_assign = [(ph, bh) for ph, pt, bh, bt in strata]
    obs = stat(obs_assign)
    ge = 0
    pool = []
    for ph, pt, bh, bt in strata:
        pool.append((pt, bt, ph + bh))  # (n_patho, n_benign, n_hits)
    for _ in range(n):
        assign = []
        for (pt, bt, hits), (ph, _pt, bh, _bt) in zip(pool, strata):
            # draw how many of `hits` land on pathogenic (hypergeometric)
            items = [1] * pt + [0] * bt
            rng.shuffle(items)
            a = sum(items[:hits])
            assign.append((a, hits - a))
        if stat(assign) >= obs - 1e-12:
            ge += 1
    return obs, (ge + 1) / (n + 1)


def main():
    strata = load()
    tot_ph = sum(s[0] for s in strata); tot_pt = sum(s[1] for s in strata)
    tot_bh = sum(s[2] for s in strata); tot_bt = sum(s[3] for s in strata)
    mh_or, chi, p_one, used = mantel_haenszel(strata)
    obs_diff, p_perm = perm_test(strata)
    from scipy.stats import fisher_exact
    f_or, f_p = fisher_exact([[tot_ph, tot_pt - tot_ph],
                              [tot_bh, tot_bt - tot_bh]], alternative="greater")
    out = Path("results/explore/variant_cmh.md")
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        f.write("# Variant-at-interface, cluster-aware re-analysis\n\n")
        f.write("Variants are stratified by complex so that within-protein "
                "correlation does not inflate significance.\n\n")
        f.write("- Complexes contributing 2x2 information: {}\n".format(used))
        f.write("- Pooled counts: pathogenic {}/{} ({:.1%}) vs benign {}/{} "
                "({:.1%}) at interface\n".format(
                    tot_ph, tot_pt, tot_ph / tot_pt, tot_bh, tot_bt, tot_bh / tot_bt))
        f.write("- Naive (pooled) Fisher, ignoring clustering: OR {:.2f}, "
                "one-sided p = {:.4f}\n".format(f_or, f_p))
        f.write("- **Mantel-Haenszel (stratified by complex): common OR {:.2f}, "
                "CMH chi2 = {:.2f}, one-sided p = {:.4f}**\n".format(
                    mh_or, chi, p_one))
        f.write("- Within-complex permutation (20,000): observed rate "
                "difference {:.3f}, one-sided p = {:.4f}\n\n".format(
                    obs_diff, p_perm))
        surv = "survives" if p_one < 0.05 else "does NOT survive"
        f.write("The pathogenic-vs-benign interface enrichment {} clustering "
                "control (CMH p = {:.4f}).\n".format(surv, p_one))
    print("MH common OR {:.2f}, CMH one-sided p = {:.4f}".format(mh_or, p_one))
    print("permutation one-sided p = {:.4f}".format(p_perm))
    print("naive Fisher OR {:.2f}, p = {:.4f}".format(f_or, f_p))
    print("Wrote", out)


if __name__ == "__main__":
    main()
