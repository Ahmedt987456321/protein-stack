# Response to peer review

Working log of actions taken against the reviewer's points. Honest record:
what was done, what it showed, what remains.

## 1. TAS/IC are not experimental evidence codes (BIGGEST ISSUE)

**Verified true.** `config.yaml` used
`evidence_codes: [EXP, IDA, IPI, IMP, IGI, IEP, HTP, HDA, HMP, HGI, HEP, TAS, IC]`.
TAS (Traceable Author Statement) and IC (Inferred by Curator) are not
experimental per the GO Consortium.

**Impact measured** (all six species GAFs, config evidence set = 988,534 rows):
- strict EXP-family: 855,277 (86.5%)
- TAS + IC: 133,257 (13.5%) -- TAS 129,475 (13.1%), IC 3,782 (0.4%)

So 13.5% of ground-truth rows are non-experimental. Not negligible; a strict-EXP
rerun is required and will change protein counts, the annotated/dark partition,
temporal cohorts, and Fmax. Action: rerun the pipeline with strict EXP-family
codes (see below), report whether every result survives.

## 9. Variant-interface statistics ignore within-protein clustering

**Addressed; result survives and strengthens.** Added scripts/45_variant_cmh.py:
a Cochran-Mantel-Haenszel test stratified by complex, plus a within-complex
permutation, so per-protein correlation cannot inflate significance.

- Naive pooled Fisher (as in the paper): OR 1.99, one-sided p = 0.0103
- Mantel-Haenszel (stratified by complex): common OR 2.76, CMH one-sided p = 0.0081
- Within-complex permutation (20,000): one-sided p = 0.0069

The pathogenic-vs-benign interface enrichment holds under clustering control, and
is in fact stronger once per-complex confounding is removed. The paper's Fisher
number will be replaced/supplemented by the CMH result.

## 3. Temporal ablation: how much does structure add prospectively?

**Done, and it materially sharpens the story (honestly).** Added a sequence-only
(SEQ) arm to scripts/19 alongside SS (seq+struct) and FULL, all leakage-
controlled, scored on the same proteins.

Prospective (temporal) top-1 precision, SEQ vs SEQ+structure, same proteins:

| branch | n paired | SEQ | SEQ+struct | delta |
|---|---:|---:|---:|---:|
| molecular function | 1667 | 0.439 | 0.448 | +0.010 |
| biological process | 1694 | 0.200 | 0.211 | +0.011 |
| cellular component | 1853 | 0.274 | 0.281 | +0.007 |

On proteins sequence already covers, structure adds only ~+0.01 prospectively.
But structure's real prospective value is COVERAGE: seq+struct predicts 2,124 MF
proteins vs 1,667 for sequence alone -- structure reaches 457 more (+27%), the
no-homolog proteins, at ~0.32 precision (MF; sequence scores 0 on them by
construction), all above the permutation null.

Honest conclusion: prospectively, structure extends coverage to sequence-less
proteins rather than lifting precision where sequence works. This is consistent
with the cross-sectional twilight-zone result and argues for reframing the story
around "structure helps where sequence fails (coverage + twilight zone)" rather
than a broad "structure-based" precision claim -- supporting the reviewer's
suggested retitle.

## Wording / factual / reference fixes (points 4, 6, 10, 11, 13)

Applied to the paper source (paper/latex/main.tex, references.bib):
- "four domains of life" -> "two domains of life" everywhere (6 species span
  Bacteria + Eukaryota only). [point 6]
- "pre-registered" -> "pre-specified" (no timestamped preregistration). [point 10]
- "unannotated" -> "experimentally unannotated"; application-set wording made
  precise. [point 11]
- "Cross-kingdom replication" -> "Cross-taxon scale-up"; limitations now state the
  scale-up reuses three original species, so it is expanded validation, not an
  independent replication. [reviewer's replication caveat]
- Evidence-code list in Methods now strict experimental (TAS/IC removed), matching
  the strict-EXP rerun. [point 1]
- Novelty reframed: intro now cites DeepFRI (Gligorijevic 2021, Nat Commun 12:3168)
  and PANDA-3D (Zhao 2024, NARGB 6:lqae094) as established structure-based
  prediction; our contribution is where structure helps/hurts + gating +
  leakage-aware evaluation, not the discovery. [point 4]
- Limitations now state the fusion is not compared against modern learned
  predictors (DeepFRI/PANDA-3D). [point 5, honest disclosure]
- LAFA reference updated from arXiv to the published version (Bioinformatics
  Advances 2026, vbag221). [point 13]
- All references validated against the published record (11 checked, all exact).

## Strict-EXP rerun (point 1) -- IN PROGRESS

config.yaml evidence_codes set to strict experimental (TAS/IC removed). Full
six-species pipeline relaunched (run_all.py). Under strict-EXP, 60,115 proteins
have >=1 experimental annotation (vs 61,576 with TAS/IC). Comparing to the backed-up
lax results (results_lax_backup/) when it completes.

## 2. Seq+Domain vs Seq+Domain+Struct ablation -- BUILT INTO THE RERUN

The strongest single stream is domains (0.575), not structure (0.437), so the
reviewer asked whether structure adds value after domains are known. Added a
Seq+Domain arm (SD, structure off) to scripts/09 and the C-minus-SD ablation
column + macro line to scripts/10. The running strict-EXP rerun computes it, so
the same rerun answers both point 1 and point 2. Result reported after the run.

## Status of remaining review points

- Point 5 (modern baselines DeepFRI/PANDA-3D): honestly disclosed as a limitation
  in the paper. Running them requires installing and executing those model stacks;
  not attempted in this pass. Framing shifted to interpretability + leakage
  control + provenance, per the reviewer's own suggestion.
- Point 7 (cluster/family-held-out split): a robustness rerun with a
  cluster-based split. Planned as a follow-up run after the strict-EXP rerun
  completes.
- Point 8 (STRING channel sensitivity): analysis noted; STRING is already
  excluded from the temporal (leakage-controlled) arm, which is the part that
  matters most.
- Point 12 (move Section 4 applications to supplementary): a structural decision
  for the author. Flagged, not done unilaterally.
- Point 14 (retitle around gating + leakage): supported by the temporal-ablation
  finding (structure's prospective value is coverage, not precision). Flagged for
  the author.

## Paper edits applied (paper/latex/, off-repo)

- Variant section now reports the CMH cluster-aware statistic (OR 2.76, p=0.008).
- Figure 3 annotation updated to the CMH result.
- All wording/reference fixes above.
- Number-dependent updates (Fmax, temporal, leakage, counts) HELD until the
  strict-EXP rerun completes, then updated consistently.

## STRICT-EXP RERUN COMPLETE -- results (the pivotal outcome)

The full six-species pipeline reran with strict experimental codes and EVERY gate
passed. Comparing strict-EXP vs the lax (TAS/IC-included) baseline:

### What survives (robust)
- Twilight-zone gain (seq vs seq+struct, MF <30%): +0.101 strict vs +0.086 lax. HOLDS.
- Fusion beats every single stream (G3): PASS. HOLDS.
- Temporal top-1 precision (leakage-controlled SS, MF): 0.425 strict vs 0.421 lax,
  permutation p=0.001. HOLDS.
- **Feature-leakage measurement: +0.223/+0.369/+0.474 (MF/BP/CC) strict vs
  +0.226/+0.370/+0.472 lax. ESSENTIALLY IDENTICAL -- the paper's most novel
  result is fully robust to strict evidence.**

### The finding that changes the thesis (reviewer point 2, confirmed)
The Seq+Domain (SD) ablation shows structure is REDUNDANT with InterPro domains:

macro-Fmax: seq 0.460, structure 0.408, domain 0.563, seq+domain 0.570,
full fusion 0.569. **Structure over seq+domain (C - SD) = -0.0008.**

Per cell, C - SD ranges only -0.008 to +0.006 across all twelve branch x bin
cells -- including the twilight zone (MF <30%: -0.0076; structure slightly HURTS
once domains are present). Domains already capture fold/family and do it better
(MF <30%: domain 0.488 vs structure 0.377).

Temporal ablation agrees: structure adds only +0.009 to +0.011 prospectively over
sequence alone, its value being coverage of no-homolog proteins, not precision.

### Honest conclusion
Structure beats SEQUENCE in the twilight zone, but does NOT beat SEQUENCE+DOMAINS
-- AlphaFold structural transfer is redundant with cheap InterPro domain
assignments for GO-term prediction. This undercuts the "structure-based" framing
and confirms the reviewer's central concern. It is, however, a useful and honest
result, and the leakage-aware evaluation + feature-leakage measurement (the
genuinely novel contribution) are fully robust. The paper should be reframed
around leakage-aware evaluation, with structure-vs-domain redundancy reported as
an honest finding rather than buried.

## PAPER REFRAMED (single six-species strict-EXP dataset)

Title: "Leakage-aware evaluation of protein function prediction, and the limits
of structural transfer" (author-approved direction).

The manuscript was rebuilt around one consistent dataset: six species, two
domains of life, strict experimental evidence, 41,157 annotated proteins. This
also dissolves the reviewer's "scale-up reuses species, not a real replication"
objection, since there is no longer a separate scale-up claim.

Changes:
- Abstract and introduction reframed: structure beats sequence below 30% but is
  redundant with InterPro domains; the leakage-aware evaluation is the central
  contribution.
- Tables 1-6 migrated to six-species strict-EXP numbers; the scale-up section and
  its gate table (old Table 7) removed.
- New in Table 2: the S, SD (sequence+domain), and C - SD columns, so the
  structure-redundancy result is visible cell by cell (C - SD in [-0.008, +0.006]).
- Figure 1 redrawn with four streams; the domain and fusion curves nearly
  coincide, showing structure adds nothing once domains are present.
- Temporal section (now the centrepiece) uses the well-powered cohort
  (SS MF 0.425, n=2,087) and the paired feature-leakage measurement (+0.22 to
  +0.47); calibration and outcome counts updated.
- Discussion and limitations rewritten; a new limitation states the redundancy
  result is about structural transfer (Foldseek), not learned structure models
  (DeepFRI/PANDA-3D), which are not compared against.
- Concordant illustrative case generalised (the prior accession was not in this
  run's dark sample); PANK4 contradiction verified still present in kg.db.
- Variant section keeps the CMH result (OR 2.76); percentages corrected to match.

Stale analysis outputs (naive baseline, species, families, horizons) were
regenerated against the six-species strict data (scripts 22 and 21). Two script
bugs fixed in passing: hardcoded "246" in script 22, and script 21 crashing on
E. coli (no archived GOA release) instead of skipping it.

The four application analyses (pseudo-enzyme, variant, module, dynamics) are
hypothesis-generation layers largely independent of the evidence-code change and
were not re-run under strict-EXP, except the variant statistics (CMH, updated).

LaTeX validated: 22/22 citations resolve, refs resolve, braces balanced, no
em/en dashes, tables column-consistent.
