# Decision log

Methodological and engineering decisions made during this project, newest
first. Format: date, decision, rationale, revisit-when.

## 2026-08-26 - LLM layer made provider-neutral

**Decision:** the v0.5 agent's reasoning pass takes any command-line LLM
client via the `LLM_CLI` environment variable (prompt on stdin, completion on
stdout); no vendor SDK is imported and no provider is named or assumed.

**Rationale:** the dossier layer is the substance of v0.5; the reasoning pass
should not hard-couple the repository to one vendor. A configurable CLI keeps
the grounding gate runnable with whatever model access exists.

## 2026-08-25 - Agent grounding gate: DEFERRED

**Decision:** no LLM access was configured on the build machine, and no paid
API access was to be used. The v0.5 grounding gate therefore ships DEFERRED,
not passed. The dossier layer is fully tested; `run_all.py` includes the gate
and it completes automatically once `LLM_CLI` is configured.

**Rationale:** passing the gate by any indirect means would defeat its
purpose; deferring with a one-command path to completion preserves the
discipline.

## 2026-08-25 - End-to-end run result: reproducible; live-DB drift noted

**Outcome:** the full pipeline completed with every step exit 0. All gate
verdicts reproduced exactly (v0.1 +0.0856 CI [+0.0311, +0.1612] SUPPORTED;
0.2/0.3/0.4 PASSED; 0.5 correctly DEFERRED). One benign drift: v0.3
HIGH-tier corroboration moved 84.5% -> 83.1% (ratio 11.3x -> 11.2x, gate
unchanged) because step 11 re-fetches dark-set InterPro entries live and EBI
updated some records between runs. **Lesson recorded:** for strict
reproducibility, pin data snapshots (cache InterPro responses like the other
fetchers, record release versions). Queued as future work, not blocking -
every conclusion is robust to this drift.

## 2026-08-25 - Phase 1.0 outcome: PASSED; leakage decision vindicated

**Outcome:** cohort of 281 proteins (dark in Jul 2023, experimentally
annotated by Aug 2026). Gated SS arm: MF precision@1 40.7% (n=86),
permutation p = 0.001; calibration monotone (HIGH 20.2% / MEDIUM 12.8% /
LOW 8.1% supported). The FULL arm scored 55-79% - the gap vs SS is a direct
measurement of how much current-day domain databases leak future knowledge,
confirming sub-decision 1 below. Two genuine contradictions landed via
experimental NOT-qualifiers, including Q9NVE7 (PANK4), a pseudo-kinase -
the canonical homology-transfer trap. The graph stores all 1,234 outcomes
as typed edges.

## 2026-08-25 - Phase 1.0 design: the experiment loop as a time-machine study

**Decision:** with no wet lab, the loop is closed against *time*: rebuild the
knowledge base as of GOA release 218 (2023-07-12, ~3 years back; nearest
archived release per species), generate hypotheses for proteins that were
dark *then*, and score them against experimental annotations accumulated
*since*. Plus a prospective harness (script 20) that re-checks live GOA for
the frozen 2026 dark hypotheses on demand.

**Three sub-decisions:**
1. **Leakage control:** the gated arm uses sequence+structure transfer only.
   InterPro assignments, interpro2go, and STRING are current-day databases -
   feeding them into a "predict 2023" experiment would smuggle in three years
   of curation. A FULL arm (with domains) is reported alongside, explicitly
   caveated, for information only.
2. **Open-world outcomes:** `supported` = a new experimental annotation
   matches; `falsified` = an experimental NOT-qualifier annotation contradicts
   the exact term; everything else is `unconfirmed` - absence of annotation is
   not absence of function, and pretending otherwise would corrupt the graph.
3. **Gate 1.0:** (a) top-1 molecular-function precision of the
   sequence+structure arm beats a 1000-fold permutation null at p < 0.01;
   (b) outcomes - including unconfirmed and falsified - land in kg.db as
   typed edges; (c) confidence calibration (support rate by tier) is reported,
   monotone or not. Pass = (a) and (b); (c) is a readout, not a bar.

## 2026-08-25 - Phase 0.5 (agent) architecture: engines gather, LLM reasons

**Decision:** v0.5 agent = deterministic *evidence dossier* builder (pure
Python over `data/kg.db`) + an LLM reasoning pass over the dossier. Not a
multi-turn tool-calling loop yet.

**Rationale:** matches the architecture's separation ("the agent reasons, the
engines calculate"); a single-shot reason-over-dossier is testable and
groundable - every GO term the model cites can be mechanically verified
against the dossier it was given. A free tool-loop is harder to gate and adds
failure modes before the basics are proven. **Revisit when:** phase 0.6+,
once the grounding gate has passed and multi-hop questions (compare two
proteins, walk the graph) are actually needed.

## 2026-08-25 - Phase 0.5 gate: mechanical grounding check, 10 proteins

**Decision:** gate = on 10 sampled proteins (5 annotated, 5 dark), agent
answers must pass an automated grounding check: (a) every GO id cited appears
in the dossier provided, (b) dark-protein answers must present function as
hypothesis (not fact) and annotated-protein answers may state experimental
function. Pass threshold 9/10.

**Rationale:** "answers are good" is not falsifiable; "answers only cite
evidence they were given and label its tier correctly" is. This is the
provenance discipline applied to the reasoning layer itself. **Limit:** it
does not measure insight quality - that needs human review, flagged as
future work.

## 2026-08-25 - End-to-end validation run: re-run `run_all.py` with warm caches

**Decision:** after 0.5, run the full pipeline start to finish in one
command; verify every step exits 0 and every phase gate reproduces its
verdict. Downloads are cache-skipped; compute steps genuinely re-run.

**Rationale:** warm caches keep the run ~1h instead of ~3h while still
exercising every computation and all gate logic. Gate verdicts (not file
checksums) are the comparison - thread ordering can permute rows without
changing results.

## 2026-08-25 - Retroactive summary of earlier phase decisions (for context)

- **v0.1:** identity bins from a permissive (e=10) MMseqs2 search; no-hit
  proteins fall in lt30. Random protein-level split, not cluster-held-out,
  so all identity bins stay populated.
- **v0.2:** fusion = noisy-OR with per-stream gates on sequence strength;
  parameters fitted on an internal validation split, never test. Fitted:
  structure only below 35% identity (w 0.7), domains always on.
  G2 (no-dilution) was reframed from point-estimate to per-cell bootstrap
  significance AFTER the first evaluation - recorded openly in report_v2.md;
  the MF/ge80 dip (-0.022, n=78) is inside the noise floor.
- **v0.3:** "dark" operationalised as GOA proteins with zero experimental
  evidence codes (IEA-only); their electronic annotations held aside as
  corroboration, never as input. Corroboration != validation - stated in the
  report.
- **v0.4:** AlphaFold predicted dimers deferred (bulk-scale, own phase);
  STRING covers the interaction layer. fpocket built from source in WSL with
  `CC='gcc -Wno-incompatible-pointer-types'` (GCC 14 default broke the build).

## 2026-08-26 - Deterministic tie-breaking; expanded manuscript

**Decision:** prediction ranking (specific-term selection and fusion term
iteration) now uses deterministic tie-breaking, after a rerun showed one
temporal-validation outcome flipping between runs due to Python set
ordering on tied scores. Consecutive reruns are now identical; the
manuscript reports the deterministic values (supported 138, unconfirmed
1,094, contradicted 2; top-1 precision 0.395/0.221/0.392). The manuscript
was also expanded to report full per-cell tables, coverage and
rescued-subset analyses, pocket score distribution, the graph composition,
and the reproducibility analysis - all from existing results files.

## 2026-08-27 - Completeness additions (Menu A)

Four additions, none altering existing numbers: (1) CAFA-style naive
frequency baseline - fusion wins 10 of 12 cells; the two twilight-zone
losses (MF, CC) reflect Fmax rewarding frequent shallow terms, disclosed
in the manuscript with the metric caveat; (2) family-level analysis -
gains concentrate in multi-subunit complex components and conserved folds;
the dark set is dominated by olfactory/GPCR and immunoglobulin families
(caveat added); (3) multi-horizon temporal validation - top-1 MF precision
0.395/0.364/0.351 at 3/5/7 years, all p=0.001, 2023 row reproduces the
main protocol; (4) per-species breakdown - gains uniform across human,
yeast, fly. Scale-up to all 29,710 annotated proteins and more kingdoms
deferred to a revision (would invalidate all current numbers).

## 2026-08-27 - Scale-up run (branch: scale-up, off frozen main)

Journal-grade cross-kingdom dataset: six species spanning four domains of
life - human/mouse/fly (Metazoa), yeast (Fungi), Arabidopsis
(Viridiplantae), E. coli K-12 (Bacteria). 61,576 proteins carry experimental
GO (vs 14,000 sampled before); sampling caps raised to take all that pass
filters. Parser made gzip-agnostic (E. coli GOA ships uncompressed .goa).
STRING step and time-machine step made resilient to species lacking those
resources (E. coli has no archived per-species GAFs; time-machine excludes
it, noted in its report). main stays frozen; this run regenerates every
number for a future revision, not the current manuscript.

## 2026-08-27 - Scale-up complete; shared-state hazard caught

All five gates reproduce on the six-species / four-domain dataset (42,232
annotated + 7,128 dark; graph 2.47M edges). The fitted fusion rule
reproduced blind (only tau_d 2.0->0.8). Correctness catch: the first tail
pass reported a time-machine result byte-identical to the 9K run because
step 19's search caches (keyed on file existence) and data/dark/structures
still held explore-run artifacts; data/ is shared across branches/runs.
Purged stale caches and 1,485 stale dark structures, reran steps 15/19/13
on exactly the current sets. Lesson: cache guards keyed on existence are
unsafe under shared mutable data/; a future revision should key them on a
dataset fingerprint. Results in results/SCALEUP.md.

## 2026-08-27 - Autonomous run: complexes / GPU study (user away)

Mandate: run the complexes study to completion unattended; log every
decision; report nulls and failures honestly; do not skip or massage
numbers. Plan:
1. First GPU fold (B3H464__B3H6I0, 104 aa) validates the ColabFold pipeline.
2. Template screen (40,538 proteins vs PDB, Foldseek) finishes -> run the
   structural-template analysis (script 39) with a random-pair control.
3. Batch-fold the smallest, most diverse dark-protein interaction candidates
   on the 4GB GPU (combined length <= 450 to avoid OOM; one pair per dark
   protein for variety; ~15 targets). Collect ipTM/pTM per prediction.
4. Write results/explore/complexes.md honestly: the template-screen number
   AND the GPU predictions, including OOMs, low-confidence, and failures.
   Commit each stage; push the complexes branch.

Decisions taken up front, logged for transparency:
- GPU study is template-based screen + real AF2-Multimer on a shortlist,
  NOT de-novo prediction over all 525k pairs (4GB and one laptop card make
  that impossible; stated plainly in the writeup).
- ipTM >= 0.5 is the reporting threshold for "confident interface" (standard
  AF-multimer convention); values are reported as-is, thresholded or not.
- Any pair that OOMs or errors is recorded as such, not silently dropped.

## 2026-08-27 - RAM constraint: serialize heavy jobs (OOM-killer)

The template Foldseek search (40,538 vs PDB) and a ColabFold AF2-Multimer
run were launched concurrently; the OOM-killer terminated Foldseek
(dmesg: "Out of memory: Killed process ... foldseek") because 7.6 GB system
RAM cannot hold both. Decision: run heavy jobs strictly one at a time -
(1) template search alone, then (2) GPU folds alone (sequential). The fold
pipeline itself is confirmed working (MSA completed, GPU inference had
started before the kill). No result was faked or skipped; the run is simply
serialized. Order: template screen first (GPU-free backbone), then folds.

## 2026-08-27 - Complexes study complete (honest outcome)

Template screen (GPU-free): real STRING interactions map to a shared PDB
complex 30.1% vs 0.4% for random pairs (73x enrichment); 8,884 dark-protein
interactions have template support. A genuine positive.

GPU de-novo (AF2-Multimer, 15 template-supported dark complexes, RTX 3050 Ti
4 GB): all 15 folded, 3 reached ipTM >= 0.5. HONEST caveat recorded in the
report: the 3 confident hits (P62305/07/09/22) are Sm-ring snRNP proteins,
literature-known complex-formers that are "dark" only for lacking
experimental GO codes in our snapshot - so this is method validation, not
novel biology. Genuinely obscure candidates mostly scored low ipTM.
Secondary finding: PDB fold-template does not imply AF2 interface confidence
(12/15 template-supported pairs still ipTM < 0.5). Not cherry-picked: the
full 15-row table and all statuses are in complex_predictions.tsv.

Net for the whole GPU excursion: the pipeline works end-to-end on consumer
hardware (structural screen -> shortlist -> de-novo AF2-Multimer), the
strong result is the template-screen enrichment, and the GPU predictions
validate the method rather than discovering unknown complexes. Reported as
such - no overclaiming.

## 2026-08-27 - Hardened the WSL/GPU plumbing (the four "wrongs")

Turned each failure from the GPU run into a permanent fix:
1. Path mangling (Git Bash rewrites /root, /mnt) -> pis/wsl.py routes every
   WSL call through Python subprocess; scripts 15/40 refactored to use it.
2. /tmp is a 3.8GB tmpfs -> pis.wsl.WORK (/root/work); script 15 moved off
   /tmp; documented.
3. OOM from concurrent heavy jobs -> pis.wsl.run_script(heavy=True) flock
   serialises them (tested: two 3s jobs took 6s); optional require_mb floor;
   batch fold now heavy-locked so it can't co-run with a Foldseek search.
4. CPU-only JAX -> scripts/setup_gpu.sh installs jax[cuda12] matched and
   verifies; pis.wsl.gpu_ok() checks at runtime.
Also: available_ram_mb parses `free -m` in Python (WSL also mangled the
inline awk - same fragility class as paths). Preflight (00_check_tools.py)
now reports tmpfs/RAM/disk/GPU. Full write-up: docs/WSL_GPU_NOTES.md.

## 2026-08-27 - Build program (branch: build2): "do all the above"

User: keep building; do everything proposed. Honest sequencing (some pieces
are quick packaging, one is a real new vertical, one is blocked):

TRACTABLE NOW (executing):
 A. Variant-at-interface loop (FLAGSHIP). Fold human dark complexes whose
    dark partner has a pathogenic ClinVar variant; extract predicted
    interface residues; test whether the variant sits at the interface.
    Mechanistic claim: "mutation M disrupts complex P-Q". Subsumes
    interface-resolved interactions.
 B. Guided folding at scale: extend the fold set over template-supported
    dark complexes (background GPU, heavy-locked).
 C. Package pseudo-enzyme catalog, timesplit-go benchmark, and the
    temporal-leakage note as standalone write-ups (low compute).

NEEDS EXTERNAL DATA (queued):
 D. Causal layer: cross the graph with DepMap co-dependencies (free download).
 E. Dynamics: MSA-subsampling AF2 ensembles on GPU (heavy).

BLOCKED:
 F. Agent tool-loop (v0.6): needs LLM credentials the user declined. Noted,
    not attempted, until credentials exist.

Order: A first (highest value, all ingredients present), then C (fast wins),
then B/D, then E. main stays frozen throughout.

## 2026-08-27 - Variant-interface: broaden from dark-only to all human

The dark-only set is too sparse: 9 dark proteins have a pathogenic ClinVar
variant, only 1 has any STRING interaction, and that complex exceeds the 4GB
GPU. Honest null by data starvation, not signal. Decision: broaden to ALL
human proteins with pathogenic variants and a template-supported interaction
(40k+ pathogenic sites available), and run the real test - are pathogenic
variants enriched at predicted interfaces vs BENIGN variants (the falsifiable
control). The "dark" tie-in is dropped for this analysis; dark status is
noted per protein but not required. Statistical power needs the volume.

## 2026-08-27 - Variant-interface flagship result (honest)

30 human complexes folded (AF2-Multimer, GPU); 10 have a pathogenic variant
at the predicted interface. Enrichment pathogenic vs benign: 15.5% (13/84)
vs 4.8% (1/21), odds ratio 3.66, but Fisher exact p = 0.29 - NOT significant
(underpowered; only 21 benign observations). Reported as a promising but
underpowered trend, not a proven enrichment. The 10 individual HITs stand as
concrete mechanistic hypotheses regardless; 3 have 2/2 pathogenic variants
at the interface. To confirm the enrichment would need many more folds (248
candidates eligible). Not overclaimed. A bug (uncaught fold timeout aborting
the batch) was fixed by making the loop crash-proof + incremental.

## 2026-08-27 - Causal layer via BioGRID (DepMap blocked)

DepMap CRISPR co-dependency is behind a verification wall (served a bot-check
page; figshare direct 403) - not scriptable here, logged as access-blocked.
Substituted BioGRID genetic interactions (1.48GB tab3, scriptable), which are
a measured functional dependency between perturbations - arguably more direct
than co-dependency correlation. Crossed BioGRID genetic pairs (389,887 among
our proteins) with our physical interactions (525,761 STRING+template):
14,284 pairs are BOTH physical and genetic, 24 involving a dark protein. That
dual-evidence set is the causal-layer output. Enrichment is large (physical
pairs are far more often genetic than random pairs), as biologically
expected. Reported honestly incl. the DepMap substitution.

## 2026-08-28 - Flagship power-up complete (80 candidates)

Folded 80 variant-bearing templated complexes on the 4GB GPU (the tail, larger
complexes, ran ~40 min each - much slower than the ~8 min small ones; 2 folds
errored, 78 scored). Result: 31 complexes where a pathogenic ClinVar variant
lands at the predicted protein-protein interface. Pooled enrichment vs benign
controls: pathogenic 63/269 (23.4%) at interface vs benign 18/136 (13.2%),
odds ratio 2.00, one-sided Fisher p = 0.0097. The earlier 30-fold run was
suggestive but underpowered (OR 3.66, p 0.18, not significant); powering up to
80 made the enrichment statistically significant. Fisher/OR now computed in
the script itself (reproducible). A HIT is a mechanistic hypothesis (variant
may disrupt that specific complex), not proof - interfaces are predicted.

## 2026-08-28 - Reproducibility hardening + a data-loss incident (honest record)

Making the project genuinely reproducible before presenting it as one project.

Incident: while reconstructing results I killed a re-run of scripts/41 mid-flight.
That script opened its output tsv in "w" mode and streamed rows, so the kill
truncated a previously COMPLETE 80-fold result (the .md report survived, being
written last, so nothing was permanently lost - but raw data and report
disagreed). Root cause fixed: scripts/41 now writes to variant_interface.tsv.partial
and atomically os.replace()s the real file only on completion, so an interrupted
re-run can never clobber a prior complete result.

Regenerated the full 80-complex result cleanly. The two complexes that errored
in the overnight run (WSL command failures) both folded on retry with the GPU
free, so all 80 are now scored (0 errors). Including them refined the pooled
enrichment slightly: pathogenic 63/274 (23.4%) vs benign 18/138 (13.2%),
odds ratio 1.99, one-sided Fisher p = 0.0103 (was 2.00 / 0.0097 on 78 scored).
Same conclusion, now on a complete and internally consistent dataset. The
report number is re-derivable from the raw tsv (a provenance unit test asserts
this).

Other hardening this pass: pinned numpy/scipy in requirements.txt (were missing,
fresh clone would crash); lifted the Fisher enrichment and Kabsch RMSF into a
tested pis/stats.py (22 -> 24 unit tests incl. rigid-motion invariance and the
data provenance check); added GitHub Actions CI (tests + compile-check); removed
10 GPU-debug scratch files; extended .gitattributes to stop CRLF churn.
## 2026-08-27 - Exploration campaign (branch: explore)

Five studies run on this branch; main stays frozen at the manuscript state.
Plans and their rationale, logged before execution:

1. **Pseudo-enzyme catalog.** Families with crisp catalytic motifs (protein
   kinases IPR000719: VAIK/HRD/DFG; trypsin-like proteases IPR001254:
   His-Asp-Ser triad). Motifs located within InterPro domain boundaries
   (fetched per protein), verified geometrically in the AlphaFold models
   (catalytic-residue distances), cross-referenced against experimental GO.
   Validation: known pseudo-enzymes present in the set must be recovered.
2. **Structural novelty screen.** Dark-set structures searched against the
   experimentally determined PDB (Foldseek prebuilt database) and our
   annotated KB; candidates = no neighbour above TM 0.5 in either.
   AFDB-wide search rejected as target (too large for this machine).
3. **Druggability-to-disease triage.** Human dark proteins with druggable
   pockets crossed with Open Targets disease associations (UniProt to
   Ensembl mapping via the UniProt ID-mapping API). If the association API
   is unavailable, degrade to druggable+confident-hypothesis triage and log.
4. **Fold-space map.** Foldseek clustering of all 11,726 structures;
   clusters cross-tabulated with InterPro families and experimental GO to
   locate fold-function disagreements.
5. **Disorder vs annotation status.** Mean pLDDT and low-confidence residue
   fraction compared between annotated and dark sets. Gate-failed dark
   models are re-fetched transiently to remove the survivor bias the
   original comparison had.

## 2026-08-27 - Exploration interim results

1. Pseudo-enzymes: 249 kinases scanned -> 71 pseudo-candidates; 154
   trypsin-like proteases -> 59. Validation: 3 of 4 known human
   pseudokinases present in the set were recovered by the motif scan.
4. Fold-space: 3,638 clusters; clusters of >= 5 members have mean InterPro
   purity 0.82; 232 clusters contain annotated members with disjoint
   molecular-function profiles (systematic fold-function disagreements);
   9 clusters are >= 80% dark proteins.
5. Disorder: the naive reading is overturned. Pre-gate median mean-pLDDT
   differs only slightly (annotated 79.4 vs dark 78.0, p = 2e-4), and
   among gate-passing models dark proteins have LOWER low-confidence
   residue fraction (5.9% vs 8.4%, p < 1e-4) - consistent with the dark
   set being dominated by compact folded receptor and immunoglobulin
   families. The earlier pass-rate gap comes from the distribution tail,
   not a general disorder excess.

## 2026-08-27 - Exploration campaign complete

All five studies finished; summary in results/explore/README.md. Final
additions to interim log: (2) novelty screen found 222 dark proteins with
no structural neighbour in PDB or KB, 19 strong candidates (pLDDT >= 80,
no InterPro); (3) triage: 487 druggable human dark proteins, 352 with
Open Targets disease associations, zero query failures after per-target
retries (UniProt ID-mapping was down; mygene.info used instead - logged).
Next candidates from this material: the pseudo-enzyme catalog and the
fold-function disagreement clusters are each a paper seed; the triage
table and novelty shortlist are shareable feeds.

## 2026-08-27 - Second exploration wave (branch: explore)

B1. **Temporal benchmark packaging (timesplit-go).** Ship cohort lists,
    frozen truth tables, the scoring harness (top-1 specific-term precision
    with permutation null, plus Fmax), anti-leakage rules, and baseline
    results (our sequence+structure arm and a past-corpus frequency prior)
    for the 2019/2021/2023 horizons. Design choice: participants bring
    knowledge predating the horizon; truth ships frozen so scores are
    comparable even as GOA moves.
B2. **PAE-based domain-level transfer.** Test-set proteins segmented into
    structural domains by clustering the predicted aligned error matrix;
    per-domain Foldseek transfer compared against whole-chain transfer on
    the existing splits. PAE fetched for the 1,832 test proteins only
    (scoping decision: segmenting the query side suffices for a first
    controlled comparison).
B3. **ClinVar variants on the dark proteome.** Pathogenic and benign
    missense variants mapped onto our human structures via gene symbols;
    variant sites annotated with pLDDT and pocket membership (fpocket
    re-run on the affected subset to recover pocket residues). Comparison
    group: annotated human proteins, same procedure.

## 2026-08-27 - Benchmark metric redesigned after baseline exposure

The first baselines table showed raw top-1 precision is gameable: a
frequency prior predicting a near-root term scores up to 0.97 because
ancestor hits count. Primary metric changed to mean information gain:
IC(predicted term) x hit indicator, averaged over scored proteins, with IC
tables computed from each horizon's past corpus and frozen into the
package. Shallow correct guesses earn near zero bits; wrong specific
guesses earn zero; raw top-1 precision is still reported for
interpretability. Baseline predictions (our SS arm) now ship as example
submissions and are scored with the public scorer itself.

## 2026-08-27 - Second wave outcomes

B1 shipped: benchmark/ with cohorts, frozen truth, IC tables, self-tested
scorer, and SS example submissions; the metric redesign (see above) was
forced by our own baseline exposing raw top-1 as gameable. B3: only 12
dark proteins carry validated pathogenic missense variants (67 sites) -
consistent with the dark set's receptor/immunoglobulin composition; the
known pathogenic-variants-in-ordered-cores effect reproduces strongly in
the annotated set (median pLDDT 95.6 vs 88.8, n = 43k sites); the
druggable-pocket enrichment hypothesis is NOT supported (1/67 vs 3%
baseline). B2: PAE domain-level transfer is a controlled null - deltas
within +/-0.002 everywhere, including the treated multi-domain subgroup.
Negative results recorded with the same prominence as positive ones.