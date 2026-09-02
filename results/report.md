# v0.1 twilight-zone experiment - report

**Question:** does structure (AlphaFold + Foldseek) improve GO-term prediction where sequence similarity fails (<30% identity)?

**Primary endpoint** (molecular_function, lt30 bin, n=243):

- delta Fmax (Arm B - Arm A): **+0.0801**  (95% bootstrap CI [+0.0343, +0.1335])
- Verdict: **SUPPORTED**

## Fmax by arm / branch / identity bin

| branch | bin | n | Arm A (seq) | Arm B (+struct) | delta |
|---|---|---:|---:|---:|---:|
| molecular_function | lt30 | 243 | 0.2970 | 0.3771 | +0.0801 |
| molecular_function | 30to50 | 1786 | 0.4925 | 0.4361 | -0.0564 |
| molecular_function | 50to80 | 2091 | 0.6130 | 0.4773 | -0.1356 |
| molecular_function | ge80 | 2190 | 0.6764 | 0.5865 | -0.0899 |
| biological_process | lt30 | 286 | 0.2371 | 0.2557 | +0.0187 |
| biological_process | 30to50 | 1819 | 0.3545 | 0.3036 | -0.0509 |
| biological_process | 50to80 | 2087 | 0.4456 | 0.3398 | -0.1058 |
| biological_process | ge80 | 2091 | 0.4754 | 0.3999 | -0.0755 |
| cellular_component | lt30 | 301 | 0.3353 | 0.4057 | +0.0704 |
| cellular_component | 30to50 | 1952 | 0.5043 | 0.4479 | -0.0564 |
| cellular_component | 50to80 | 2250 | 0.5994 | 0.4632 | -0.1362 |
| cellular_component | ge80 | 2334 | 0.6601 | 0.5763 | -0.0838 |

Evidence type of every prediction here: COMPUTATIONAL. Ground truth: EXPERIMENTAL/CURATED GO annotations only.
