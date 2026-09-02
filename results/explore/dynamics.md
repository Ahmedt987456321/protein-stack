# Conformational dynamics via MSA-subsampling ensembles

AlphaFold gives one static model per protein. Running it at shallow MSA depth (16:32) across 8 random seeds makes it sample alternative conformations (Del Alamo et al., eLife 2022); the spread of CA positions across the ensemble is a proxy for intrinsic flexibility that the single AlphaFold DB model hides.

Per-residue RMSF is measured on centroid-superposed CA traces over residues present in every model. Dark (unannotated) proteins were prioritised - flexibility is a functional clue for a protein with no known function.

| protein | length | models | mean RMSF (A) | max RMSF (A) | call |
|---|---|---|---|---|---|
| A8DY46 | 50 | 8 | 2.0 | 6.68 | rigid |
| A0A0B4K7B8 | 50 | 8 | 1.08 | 2.78 | rigid |
| M9NEV5 | 50 | 8 | 1.04 | 2.61 | rigid |
| A8DYL9 | 50 | 8 | 0.97 | 2.79 | rigid |
| Q8TGS2 | 50 | 8 | 0.76 | 1.84 | rigid |
| B3H6I0 | 50 | 8 | 0.36 | 0.96 | rigid |

High mean RMSF flags a candidate conformationally flexible or multi-state protein; localized high-RMSF stretches (see dynamics_rmsf.tsv) are candidate hinges or disordered segments. RMSF here is rotation-free (centroid-aligned), a relative-flexibility proxy, not an absolute B-factor.
