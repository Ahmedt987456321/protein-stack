# Genetic x physical interaction crossing (causal layer)

DepMap CRISPR co-dependency is not scriptable from this environment (verification wall); BioGRID genetic interactions are used instead - a measured functional dependency between two genes' perturbations.

- Physical interaction pairs (STRING >= 700 + AF-template): 525761
- Distinct genetic-interaction pairs among our proteins (BioGRID): 389887
- Pairs that are BOTH physical and genetic: **14284** (2.7% of physical pairs also genetically interact)
- Of those, involving a dark (unannotated) protein: **24** - these carry both physical and functional-genetic evidence for a protein with no experimental annotation.

A pair with both evidence types is a stronger, more causal functional link than either alone: they physically associate AND perturbing one modifies the other's phenotype. Full list: causal_pairs.tsv.

## Dark-protein pairs with both physical and genetic evidence (top 20)

| dark-involving pair | genetic assay(s) |
|---|---|
| P53823 - Q03144 | Phenotypic Enhancement |
| P43544 - P53823 | Phenotypic Enhancement |
| P38150 - P40361 | Synthetic Growth Defect |
| P0AAA7 - P76372 | Negative Genetic |
| P27243 - P76372 | Negative Genetic |
| P37313 - P77268 | Negative Genetic |
| P07021 - P0A910 | Positive Genetic |
| P77308 - Q47622 | Positive Genetic |
| P0AAG0 - Q47622 | Negative Genetic |
| P0AFH2 - Q47622 | Negative Genetic |
| P36646 - P36678 | Negative Genetic |
| Q91VM3 - Q9CQY1 | Negative Genetic |
| Q16222 - Q3KQV9 | Negative Genetic |
| P0ACH5 - P31449 | Negative Genetic |
| P37641 - P37671 | Negative Genetic |
| P77171 - P77309 | Negative Genetic |
| P23917 - P76586 | Negative Genetic |
| P09378 - P31449 | Positive Genetic |
| P0ACI0 - P77396 | Positive Genetic |
| P0ACQ7 - P77309 | Positive Genetic |
