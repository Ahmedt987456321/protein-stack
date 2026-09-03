# v0.4 - interactions + pockets - report

STRING edges (combined score >= 700): 513431 among our 41157 proteins. Interaction arm I = guilt-by-association transfer from STRING partners; C+I = noisy-OR on top of the v0.2 fused arm with w_i = 0.2 fitted on the v0.2 validation split.

fpocket ran over 7054 dark-protein models; 4378 have a pocket with druggability >= 0.5 (see results/druggability_top.md).

## Gate

- macro-Fmax: C 0.5691 -> C+I 0.5762 (PASS); I alone 0.4006
- dilution: PASS

**Phase 0.4 gate: PASSED**

## Fmax by branch / bin

| branch | bin | n | I | C | C+I | delta |
|---|---|---:|---:|---:|---:|---:|
| molecular_function | lt30 | 226 | 0.3722 | 0.4902 | 0.5000 | +0.0098 |
| molecular_function | 30to50 | 1700 | 0.3556 | 0.6443 | 0.6489 | +0.0046 |
| molecular_function | 50to80 | 2013 | 0.3364 | 0.6746 | 0.6753 | +0.0006 |
| molecular_function | ge80 | 2184 | 0.3559 | 0.6641 | 0.6650 | +0.0009 |
| biological_process | lt30 | 237 | 0.4018 | 0.4388 | 0.4495 | +0.0107 |
| biological_process | 30to50 | 1687 | 0.3891 | 0.4775 | 0.4883 | +0.0108 |
| biological_process | 50to80 | 1954 | 0.3249 | 0.4927 | 0.4984 | +0.0056 |
| biological_process | ge80 | 1947 | 0.3198 | 0.4689 | 0.4724 | +0.0034 |
| cellular_component | lt30 | 264 | 0.5261 | 0.5083 | 0.5286 | +0.0203 |
| cellular_component | 30to50 | 1829 | 0.4953 | 0.6465 | 0.6539 | +0.0074 |
| cellular_component | 50to80 | 2063 | 0.4611 | 0.6800 | 0.6856 | +0.0055 |
| cellular_component | ge80 | 2225 | 0.4690 | 0.6433 | 0.6481 | +0.0048 |

Cells improved by > 0.005: molecular_function/lt30 (+0.0098), biological_process/lt30 (+0.0107), biological_process/30to50 (+0.0108), biological_process/50to80 (+0.0056), cellular_component/lt30 (+0.0203), cellular_component/30to50 (+0.0074), cellular_component/50to80 (+0.0055)

Evidence type of interaction transfer: COMPUTATIONAL (STRING combined scores include predicted channels).
