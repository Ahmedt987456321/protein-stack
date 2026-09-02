# PAE domain-level structural transfer

416 of 1832 test proteins are multi-domain by PAE segmentation (windows 10 residues, cut 10.0 A, min domain 40 residues); their structure stream is replaced by unioned per-domain transfer. Single-domain proteins are untouched, so all differences below come from the treated subgroup.

## all test proteins

| branch | n | C | C-dom | delta |
|---|---|---|---|---|
| molecular_function | 1514 | 0.6417 | 0.6416 | -0.0002 |
| biological_process | 1443 | 0.4912 | 0.4911 | -0.0001 |
| cellular_component | 1590 | 0.6401 | 0.6400 | -0.0001 |

## multi-domain subgroup

| branch | n | C | C-dom | delta |
|---|---|---|---|---|
| molecular_function | 355 | 0.6081 | 0.6067 | -0.0014 |
| biological_process | 326 | 0.4571 | 0.4591 | +0.0020 |
| cellular_component | 372 | 0.6196 | 0.6190 | -0.0006 |


## Interpretation

A controlled null result. Per-domain search neither helps nor hurts, even
restricted to the multi-domain subgroup: Foldseek's whole-chain alignments
already recover domain-local matches, and the fusion gate suppresses the
structure stream in the identity range where most multi-domain proteins
have sequence neighbours anyway. The hypothesis that whole-chain transfer
loses domain-local function signal is not supported at this dataset scale.
