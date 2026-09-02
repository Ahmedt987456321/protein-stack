# v1.0 - the experiment loop (time-machine study) - report

Past world: human goa_human.gaf.218.gz (2023-07-12); yeast goa_yeast.gaf.119.gz (2023-07-12); fly goa_fly.gaf.119.gz (2023-07-12); mouse goa_mouse.gaf.204.gz (2023-07-12); arabidopsis goa_arabidopsis.gaf.191.gz (2023-07-12). Cohort: **5512 proteins** that had zero experimental GO then and have it now; past KB: 36720 proteins (annotations truncated to the past release).

Arm SS = sequence+structure transfer only (leakage-controlled, gated). Arm FULL adds current-day InterPro domains - reported with that caveat, not gated.

## Top-1 precision vs permutation null (n = scored proteins)

| arm | branch | n | precision@1 | perm. p |
|---|---|---:|---:|---:|
| SS | molecular_function | 2124 | 0.421 | 0.0010 |
| SS | biological_process | 2144 | 0.188 | 0.0010 |
| SS | cellular_component | 2284 | 0.267 | 0.0010 |
| FULL | molecular_function | 3149 | 0.644 | 0.0010 |
| FULL | biological_process | 3260 | 0.592 | 0.0010 |
| FULL | cellular_component | 3140 | 0.756 | 0.0010 |

## Outcomes written back (arm SS, top 3 per branch)

supported 3538 | contradicted 10 | unconfirmed 23063 - all stored as typed edges in kg.db (rebuild step 13). `unconfirmed` is weak negative evidence only: GO is open-world.

## Confidence calibration (arm SS, supported rate by tier)

| tier | supported | total | rate |
|---|---:|---:|---:|
| HIGH | 883 | 4770 | 18.5% |
| MEDIUM | 1292 | 9898 | 13.1% |
| LOW | 1363 | 11943 | 11.4% |

Calibration monotone (HIGH >= MEDIUM >= LOW): **yes**

**Gate 1.0** (MF precision@1 of arm SS beats permutation null at p < 0.01, and outcomes stored): **PASSED**
