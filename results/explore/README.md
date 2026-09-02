# Exploration campaign - five studies on the AlphaFold dataset

Branch: explore. Main is frozen at the manuscript state; nothing here
changes any number in the paper. Scripts 30-34; per-study reports in this
directory; decisions in DECISIONS.md.

## 1. Pseudo-enzyme catalog (pseudoenzymes.md)

Motif scans within InterPro domain boundaries, geometry-checked in the
AlphaFold models. 249 protein kinases -> 71 pseudokinase candidates; 154
trypsin-like proteases -> 59 pseudo-protease candidates. Validation: 3 of
the 4 known human pseudokinases present in the dataset were recovered.
Candidates that nonetheless carry experimental catalytic GO annotations are
flagged as review-worthy conflicts.

## 2. Structural novelty screen (novelty.md)

Every dark-set model searched against the experimentally determined PDB
(Foldseek prebuilt database) and the annotated KB. 222 of 2,564 dark
proteins have no neighbour at TM >= 0.5 in either; 19 of those are
high-confidence models (mean pLDDT >= 80) with no InterPro entry - the
strongest candidate novel folds. Low-pLDDT low-TM rows are flagged as
possible modelling failures rather than novelty.

## 3. Druggability-to-disease triage (target_triage.md)

487 human dark proteins carry a druggable pocket (fpocket >= 0.5); 386
mapped to Ensembl genes; 352 have at least one Open Targets disease
association. Top of the table: AOPEP (druggability 0.97, dystonia 31),
SLC38A8 (foveal hypoplasia, association 0.81, high-confidence transporter
hypothesis), UBR7 (Li-Campeau syndrome). Each row joins pocket score,
disease link, and the pipeline's function hypothesis with its confidence
tier.

## 4. Fold-space map (foldspace.md)

Foldseek clustering of all 11,726 models: 3,638 clusters (2,061
singletons). Clusters with >= 5 members have mean InterPro purity 0.82 -
structural clustering independently recovers curated families. 232 clusters
contain annotated members with fully disjoint molecular-function profiles
(systematic fold-function disagreements, the PANK4 pattern at scale), and
9 clusters are >= 80% dark proteins.

## 5. Disorder vs annotation status (disorder.md)

With survivor bias removed (gate-failed dark models re-measured
transiently): the pre-gate confidence gap between annotated and dark
proteins is real but small (median mean-pLDDT 79.4 vs 78.0, p = 2e-4), and
among gate-passing models dark proteins have LOWER predicted-disorder
fraction (5.9% vs 8.4%, p < 1e-4). The dark set is not disorder-rich; it
is dominated by compact, folded, hard-to-assay families, and the pass-rate
gap comes from the distribution tail.

# Second wave

## 6. timesplit-go benchmark (benchmark/ at repo root)

Cohorts, frozen truth, per-horizon IC tables, a self-contained scorer, and
our sequence+structure submissions as format examples. Primary metric is
mean information gain (bits), adopted after our own baseline showed raw
top-1 precision is gameable by near-root predictions (a frequency prior
reaches 0.97 raw precision on cellular component but earns under 0.3 bits;
the SS arm earns 1.6-2.6 bits).

## 7. ClinVar variants (variants.md)

Reference-residue-validated missense variants mapped onto the structures.
Pathogenic sites concentrate in high-confidence cores in the annotated set
(median pLDDT 95.6 vs 88.8 benign, 43k sites). Only 12 dark proteins carry
pathogenic missense variants (67 sites), and druggable-pocket enrichment
is not supported (1/67 vs 3% baseline) - a recorded negative.

## 8. PAE domain-level transfer (pae_domains.md)

Controlled null: replacing whole-chain structural transfer with unioned
per-domain transfer (PAE segmentation, 355 multi-domain test proteins)
changes Fmax by less than 0.002 in every branch and subgroup.
