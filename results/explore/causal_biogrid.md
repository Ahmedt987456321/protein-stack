# Genetic x physical interaction crossing (causal layer)

DepMap CRISPR co-dependency is not scriptable from this environment (verification wall); BioGRID genetic interactions are used instead - a measured functional dependency between two genes' perturbations.

- Physical interaction pairs (STRING >= 700 + AF-template): 522438
- Distinct genetic-interaction pairs among our proteins (BioGRID): 391754
- Pairs that are BOTH physical and genetic: **14285** (2.7% of physical pairs also genetically interact)
- Of those, involving a dark (unannotated) protein: **33** - these carry both physical and functional-genetic evidence for a protein with no experimental annotation.

A pair with both evidence types is a stronger, more causal functional link than either alone: they physically associate AND perturbing one modifies the other's phenotype. Full list: causal_pairs.tsv.

## Dark-protein pairs with both physical and genetic evidence (top 20)

| dark-involving pair | genetic assay(s) |
|---|---|
| P38150 - P40361 | Synthetic Growth Defect |
| P31697 - P45420 | Negative Genetic |
| P37313 - P77268 | Negative Genetic |
| P45760 - Q46836 | Negative Genetic |
| P75831 - P76185 | Positive Genetic |
| P0AAG0 - P76128 | Negative Genetic |
| P0AFH6 - P76128 | Positive Genetic |
| P77308 - Q47622 | Positive Genetic |
| P0AAG0 - Q47622 | Negative Genetic |
| P0AFH2 - Q47622 | Negative Genetic |
| P45420 - Q47536 | Positive Genetic |
| P36646 - P36678 | Negative Genetic |
| Q46800 - Q46809 | Negative Genetic |
| P77337 - P77739 | Negative Genetic |
| P29704 - Q05497 | Negative Genetic |
| A2P2R3 - P43123 | Negative Genetic |
| P84180 - Q9W1U5 | Phenotypic Enhancement |
| Q9BRY0 - Q9NUM3 | Positive Genetic |
| Q91VM3 - Q9CQY1 | Negative Genetic |
| P42701 - Q8IZI9 | Positive Genetic |
