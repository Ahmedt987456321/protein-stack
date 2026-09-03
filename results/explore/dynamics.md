# Conformational dynamics via MSA-subsampling ensembles

AlphaFold gives one static model per protein. Running it at shallow MSA depth (16:32) across 8 random seeds makes it sample alternative conformations (Del Alamo et al., eLife 2022); the spread of CA positions across the ensemble is a proxy for intrinsic flexibility that the single AlphaFold DB model hides.

Per-residue RMSF is measured on Kabsch-superposed CA traces (rotation + translation removed) over residues present in every model. Dark (unannotated) proteins were prioritised - flexibility is a functional clue for a protein with no known function.

| protein | length | models | mean RMSF (A) | max RMSF (A) | call |
|---|---|---|---|---|---|
| F4J9Y1 | 50 | 8 | 4.23 | 9.44 | flexible |
| A0A0B4K6L9 | 53 | 8 | 2.16 | 7.54 | rigid |
| A0A1I9LLC9 | 51 | 8 | 1.77 | 4.1 | rigid |
| Q3E787 | 52 | 8 | 1.26 | 4.36 | rigid |
| A0A286YFK9 | 51 | 8 | 0.99 | 3.59 | rigid |
| A0A5F8MPQ1 | 52 | 8 | 0.92 | 2.17 | rigid |
| P0DQW1 | 50 | 8 | 0.29 | 1.21 | rigid |
| P0C5L4 | 52 | 8 | 0.13 | 0.55 | rigid |

High mean RMSF flags a candidate conformationally flexible or multi-state protein; localized high-RMSF stretches (see dynamics_rmsf.tsv) are candidate hinges or disordered segments. RMSF here is a Kabsch-superposed relative-flexibility proxy, not an absolute B-factor.
