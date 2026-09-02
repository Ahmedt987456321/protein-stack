# Scale-up: cross-kingdom replication

Branch: scale-up (off frozen main; the manuscript numbers are unchanged).
The full pipeline was rerun on a six-species dataset spanning four domains
of life - human, mouse, fly (Metazoa), yeast (Fungi), Arabidopsis
(Viridiplantae), and Escherichia coli K-12 (Bacteria) - to test whether the
manuscript's findings survive scale and taxonomic breadth.

## Dataset

| | manuscript | scale-up |
|---|---|---|
| species | 3 (human, yeast, fly) | 6 (+ mouse, Arabidopsis, E. coli) |
| domains of life | 2 | 4 |
| proteins with experimental GO | 14,000 sampled | 61,576 (all) |
| final annotated set (pLDDT gate) | 9,162 | 42,232 |
| dark set | 2,564 | 7,128 |
| STRING edges (>= 700) | 59,547 | 525,761 |
| knowledge graph | 531,238 edges | 2,467,203 edges / 122,412 nodes |

## Every gate holds

| gate | manuscript | scale-up | verdict |
|---|---|---|---|
| 0.1 twilight-zone dFmax (MF, <30%) | +0.086 [+0.031,+0.161] | +0.080 [+0.034,+0.134] | SUPPORTED |
| 0.2 fusion macro-Fmax vs streams | 0.591 (A .459 S .437 D .575) | 0.577 (A .474 S .410 D .575) | PASSED |
| 0.3 dark HIGH-tier corroboration | 84.5% vs 7.5% | 91.3% vs 10.5% | PASSED |
| 0.4 fusion + interactions | 0.614 | 0.596 | PASSED |
| 1.0 temporal top-1 (MF, p) | 0.395, p=0.001, n=86 | 0.421, p=0.001, n=2124 | PASSED |

## The fitted rule reproduced blind

Grid-fit on the six-species validation split returned tau=0.35, alpha=0.7,
beta=0, gamma=1.0 - identical to the manuscript's "structure only below 35%
identity at weight 0.7, domains on". The single change: the domain gate
tau_d moved 2.0 -> 0.8, switching domain evidence off above 80% identity
where sequence homology already dominates. A minor, interpretable
refinement, not a contradiction.

## Notes

- The temporal validation is now far better powered (2,124-2,284 scored
  proteins per branch vs 86-97), and molecular-function precision rose to
  0.421. Cohorts differ, so absolute values are not directly comparable to
  the manuscript; the permutation significance and the qualitative pattern
  are what carry over.
- E. coli has no archived per-species GAFs, so it is excluded from the
  temporal cohort (five species) but present everywhere else.
- Pockets: 4,378 of 7,049 dark proteins druggable (>= 0.5).
- data/ is shared mutable state on this machine; the first pass reused stale
  time-machine search caches and stale structures from the earlier
  exploration run. Those were purged and the affected steps rerun on exactly
  the six-species sets before recording any number here (see DECISIONS.md).
