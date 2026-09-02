# v0.2 - domains + identity-aware fusion - report

**Arms:** A = sequence only; B = naive max-fusion (v0.1); S = structure only; D = domains only; C = fitted fusion (v0.2)

## Gates

- **G1** (fusion beats sequence in the twilight zone, MF/lt30, n=243): delta Fmax +0.2049, 95% CI [+0.1568, +0.2605] -> PASS
- **G2** (no statistically significant dilution, per-cell bootstrap 95% CI): PASS
  - no cell had a negative point estimate
- **G3** (fusion beats every single stream, macro-Fmax): PASS
  - macro-Fmax: A 0.4742, B 0.4224, S 0.4100, D 0.5750, C 0.5766

**Phase 0.2 gate: PASSED**

## Fmax by branch / bin

| branch | bin | n | A seq | B naive | S struct | D domain | C fused | C-A |
|---|---|---:|---:|---:|---:|---:|---:|---:|
| molecular_function | lt30 | 243 | 0.2970 | 0.3771 | 0.3771 | 0.5103 | 0.5019 | +0.2049 |
| molecular_function | 30to50 | 1786 | 0.4925 | 0.4361 | 0.4361 | 0.6375 | 0.6365 | +0.1440 |
| molecular_function | 50to80 | 2091 | 0.6130 | 0.4773 | 0.4773 | 0.6765 | 0.6876 | +0.0746 |
| molecular_function | ge80 | 2190 | 0.6764 | 0.5865 | 0.5303 | 0.6785 | 0.6764 | +0.0000 |
| biological_process | lt30 | 286 | 0.2371 | 0.2557 | 0.2491 | 0.4106 | 0.4134 | +0.1764 |
| biological_process | 30to50 | 1819 | 0.3545 | 0.3036 | 0.3036 | 0.4821 | 0.4882 | +0.1337 |
| biological_process | 50to80 | 2087 | 0.4456 | 0.3398 | 0.3398 | 0.4811 | 0.5181 | +0.0725 |
| biological_process | ge80 | 2091 | 0.4754 | 0.3999 | 0.3689 | 0.4918 | 0.4754 | +0.0000 |
| cellular_component | lt30 | 301 | 0.3353 | 0.4057 | 0.3961 | 0.5323 | 0.5294 | +0.1941 |
| cellular_component | 30to50 | 1952 | 0.5043 | 0.4479 | 0.4479 | 0.6426 | 0.6430 | +0.1387 |
| cellular_component | 50to80 | 2250 | 0.5994 | 0.4632 | 0.4632 | 0.6760 | 0.6897 | +0.0903 |
| cellular_component | ge80 | 2334 | 0.6601 | 0.5763 | 0.5310 | 0.6800 | 0.6601 | +0.0000 |

Evidence types: seq/structure transfer and learned domain associations are COMPUTATIONAL; interpro2go mappings are CURATED; ground truth is EXPERIMENTAL GO only.
