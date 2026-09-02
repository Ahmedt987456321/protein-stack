# v0.2 - domains + identity-aware fusion - report

**Arms:** A = sequence only; B = naive max-fusion (v0.1); S = structure only; D = domains only; SD = sequence+domain (structure off); C = fitted fusion (v0.2)

## Gates

- **G1** (fusion beats sequence in the twilight zone, MF/lt30, n=226): delta Fmax +0.2140, 95% CI [+0.1736, +0.2780] -> PASS
- **G2** (no statistically significant dilution, per-cell bootstrap 95% CI): PASS
  - no cell had a negative point estimate
- **G3** (fusion beats every single stream, macro-Fmax): PASS
  - macro-Fmax: A 0.4599, B 0.4205, S 0.4079, D 0.5627, SD 0.5699, C 0.5691
- **incremental value of structure over sequence+domain (macro C minus SD): -0.0008** (review point 2: does structure add value after domains are known?)

**Phase 0.2 gate: PASSED**

## Fmax by branch / bin

| branch | bin | n | A seq | B naive | S struct | D domain | SD seq+dom | C fused | C-A | C-SD |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| molecular_function | lt30 | 226 | 0.2762 | 0.3768 | 0.3768 | 0.4880 | 0.4978 | 0.4902 | +0.2140 | -0.0076 |
| molecular_function | 30to50 | 1700 | 0.4848 | 0.4439 | 0.4439 | 0.6512 | 0.6523 | 0.6443 | +0.1594 | -0.0080 |
| molecular_function | 50to80 | 2013 | 0.6090 | 0.4805 | 0.4805 | 0.6622 | 0.6746 | 0.6746 | +0.0657 | +0.0000 |
| molecular_function | ge80 | 2184 | 0.6641 | 0.5883 | 0.5282 | 0.6667 | 0.6641 | 0.6641 | +0.0000 | +0.0000 |
| biological_process | lt30 | 237 | 0.2213 | 0.2639 | 0.2639 | 0.4295 | 0.4355 | 0.4388 | +0.2175 | +0.0034 |
| biological_process | 30to50 | 1687 | 0.3350 | 0.2906 | 0.2906 | 0.4687 | 0.4782 | 0.4775 | +0.1425 | -0.0007 |
| biological_process | 50to80 | 1954 | 0.4251 | 0.3225 | 0.3225 | 0.4461 | 0.4927 | 0.4927 | +0.0676 | +0.0000 |
| biological_process | ge80 | 1947 | 0.4689 | 0.3977 | 0.3600 | 0.4594 | 0.4689 | 0.4689 | +0.0000 | +0.0000 |
| cellular_component | lt30 | 264 | 0.2949 | 0.4055 | 0.4012 | 0.5050 | 0.5026 | 0.5083 | +0.2134 | +0.0057 |
| cellular_component | 30to50 | 1829 | 0.5039 | 0.4422 | 0.4422 | 0.6511 | 0.6484 | 0.6465 | +0.1425 | -0.0019 |
| cellular_component | 50to80 | 2063 | 0.5919 | 0.4700 | 0.4700 | 0.6624 | 0.6800 | 0.6800 | +0.0881 | +0.0000 |
| cellular_component | ge80 | 2225 | 0.6433 | 0.5639 | 0.5154 | 0.6622 | 0.6433 | 0.6433 | +0.0000 | +0.0000 |

Evidence types: seq/structure transfer and learned domain associations are COMPUTATIONAL; interpro2go mappings are CURATED; ground truth is EXPERIMENTAL GO only.
