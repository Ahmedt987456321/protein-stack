# Variant-at-interface, cluster-aware re-analysis

Variants are stratified by complex so that within-protein correlation does not inflate significance.

- Complexes contributing 2x2 information: 26
- Pooled counts: pathogenic 63/274 (23.0%) vs benign 18/138 (13.0%) at interface
- Naive (pooled) Fisher, ignoring clustering: OR 1.99, one-sided p = 0.0103
- **Mantel-Haenszel (stratified by complex): common OR 2.76, CMH chi2 = 5.78, one-sided p = 0.0081**
- Within-complex permutation (20,000): observed rate difference 0.099, one-sided p = 0.0069

The pathogenic-vs-benign interface enrichment survives clustering control (CMH p = 0.0081).
