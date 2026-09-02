# v0.1 twilight-zone experiment - report

**Question:** does structure (AlphaFold + Foldseek) improve GO-term prediction where sequence similarity fails (<30% identity)?

**Primary endpoint** (molecular_function, lt30 bin, n=226):

- delta Fmax (Arm B - Arm A): **+0.1006**  (95% bootstrap CI [+0.0599, +0.1552])
- Verdict: **SUPPORTED**

## Fmax by arm / branch / identity bin

| branch | bin | n | Arm A (seq) | Arm B (+struct) | delta |
|---|---|---:|---:|---:|---:|
| molecular_function | lt30 | 226 | 0.2762 | 0.3768 | +0.1006 |
| molecular_function | 30to50 | 1700 | 0.4848 | 0.4439 | -0.0409 |
| molecular_function | 50to80 | 2013 | 0.6090 | 0.4805 | -0.1285 |
| molecular_function | ge80 | 2184 | 0.6641 | 0.5883 | -0.0758 |
| biological_process | lt30 | 237 | 0.2213 | 0.2639 | +0.0426 |
| biological_process | 30to50 | 1687 | 0.3350 | 0.2906 | -0.0443 |
| biological_process | 50to80 | 1954 | 0.4251 | 0.3225 | -0.1026 |
| biological_process | ge80 | 1947 | 0.4689 | 0.3977 | -0.0712 |
| cellular_component | lt30 | 264 | 0.2949 | 0.4055 | +0.1106 |
| cellular_component | 30to50 | 1829 | 0.5039 | 0.4422 | -0.0618 |
| cellular_component | 50to80 | 2063 | 0.5919 | 0.4700 | -0.1219 |
| cellular_component | ge80 | 2225 | 0.6433 | 0.5639 | -0.0794 |

Evidence type of every prediction here: COMPUTATIONAL. Ground truth: EXPERIMENTAL/CURATED GO annotations only.
