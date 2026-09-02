# v0.4 - interactions + pockets - report

STRING edges (combined score >= 700): 525761 among our 42232 proteins. Interaction arm I = guilt-by-association transfer from STRING partners; C+I = noisy-OR on top of the v0.2 fused arm with w_i = 0.6 fitted on the v0.2 validation split.

fpocket ran over 8537 dark-protein models; 5301 have a pocket with druggability >= 0.5 (see results/druggability_top.md).

## Gate

- macro-Fmax: C 0.5766 -> C+I 0.5962 (PASS); I alone 0.4009
- dilution: PASS
  - molecular_function/50to80: point -0.0019, CI [-0.0071, +0.0028] -> not significant

**Phase 0.4 gate: PASSED**

## Fmax by branch / bin

| branch | bin | n | I | C | C+I | delta |
|---|---|---:|---:|---:|---:|---:|
| molecular_function | lt30 | 243 | 0.3689 | 0.5019 | 0.5424 | +0.0405 |
| molecular_function | 30to50 | 1786 | 0.3394 | 0.6365 | 0.6409 | +0.0044 |
| molecular_function | 50to80 | 2091 | 0.3311 | 0.6876 | 0.6857 | -0.0019 |
| molecular_function | ge80 | 2190 | 0.3406 | 0.6764 | 0.6826 | +0.0062 |
| biological_process | lt30 | 286 | 0.4153 | 0.4134 | 0.4496 | +0.0362 |
| biological_process | 30to50 | 1819 | 0.3930 | 0.4882 | 0.5175 | +0.0293 |
| biological_process | 50to80 | 2087 | 0.3369 | 0.5181 | 0.5212 | +0.0031 |
| biological_process | ge80 | 2091 | 0.3224 | 0.4754 | 0.4909 | +0.0155 |
| cellular_component | lt30 | 301 | 0.5275 | 0.5294 | 0.5925 | +0.0631 |
| cellular_component | 30to50 | 1952 | 0.4985 | 0.6430 | 0.6650 | +0.0219 |
| cellular_component | 50to80 | 2250 | 0.4611 | 0.6897 | 0.6983 | +0.0086 |
| cellular_component | ge80 | 2334 | 0.4762 | 0.6601 | 0.6676 | +0.0075 |

Cells improved by > 0.005: molecular_function/lt30 (+0.0405), molecular_function/ge80 (+0.0062), biological_process/lt30 (+0.0362), biological_process/30to50 (+0.0293), biological_process/ge80 (+0.0155), cellular_component/lt30 (+0.0631), cellular_component/30to50 (+0.0219), cellular_component/50to80 (+0.0086), cellular_component/ge80 (+0.0075)

Evidence type of interaction transfer: COMPUTATIONAL (STRING combined scores include predicted channels).
