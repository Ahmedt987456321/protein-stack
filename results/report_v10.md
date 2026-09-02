# v1.0 - the experiment loop (time-machine study) - report

Past world: human goa_human.gaf.218.gz (2023-07-12); yeast goa_yeast.gaf.119.gz (2023-07-12); fly goa_fly.gaf.119.gz (2023-07-12); mouse goa_mouse.gaf.204.gz (2023-07-12); arabidopsis goa_arabidopsis.gaf.191.gz (2023-07-12). Cohort: **5462 proteins** that had zero experimental GO then and have it now; past KB: 35695 proteins (annotations truncated to the past release).

Arm SS = sequence+structure transfer only (leakage-controlled, gated). Arm FULL adds current-day InterPro domains - reported with that caveat, not gated.

## Top-1 precision vs permutation null (n = scored proteins)

| arm | branch | n | precision@1 | perm. p |
|---|---|---:|---:|---:|
| SS | molecular_function | 2087 | 0.425 | 0.0010 |
| SS | biological_process | 1964 | 0.170 | 0.0010 |
| SS | cellular_component | 2109 | 0.254 | 0.0010 |
| FULL | molecular_function | 3271 | 0.643 | 0.0010 |
| FULL | biological_process | 3150 | 0.578 | 0.0010 |
| FULL | cellular_component | 3077 | 0.751 | 0.0010 |

## Leakage, paired (SS vs FULL on the same proteins)

Both arms scored on the proteins they both predict, so the difference is leakage from present-day curated evidence, not a coverage artefact.

| branch | n paired | SS prec@1 | FULL prec@1 | leakage |
|---|---:|---:|---:|---:|
| molecular_function | 2087 | 0.425 | 0.648 | +0.223 |
| biological_process | 1964 | 0.170 | 0.538 | +0.369 |
| cellular_component | 2109 | 0.254 | 0.728 | +0.474 |

## Prospective contribution of structure (SEQ vs SS, same proteins)

Leakage-controlled: does structure add prospective signal over sequence alone?

| branch | n paired | SEQ prec@1 | SEQ+struct prec@1 | delta |
|---|---:|---:|---:|---:|
| molecular_function | 1625 | 0.440 | 0.451 | +0.011 |
| biological_process | 1508 | 0.183 | 0.193 | +0.010 |
| cellular_component | 1681 | 0.264 | 0.272 | +0.009 |

## Outcomes written back (arm SS, top 3 per branch)

supported 3280 | contradicted 11 | unconfirmed 21112 - all stored as typed edges in kg.db (rebuild step 13). `unconfirmed` is weak negative evidence only: GO is open-world.

## Confidence calibration (arm SS, supported rate by tier)

| tier | supported | total | rate |
|---|---:|---:|---:|
| HIGH | 763 | 4264 | 17.9% |
| MEDIUM | 1202 | 9057 | 13.3% |
| LOW | 1315 | 11082 | 11.9% |

Calibration monotone (HIGH >= MEDIUM >= LOW): **yes**

**Gate 1.0** (MF precision@1 of arm SS beats permutation null at p < 0.01, and outcomes stored): **PASSED**
