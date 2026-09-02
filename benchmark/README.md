# timesplit-go: a temporal benchmark for protein function prediction

Most function-prediction evaluations use a random held-out split of one
annotation snapshot. This benchmark instead asks the question practitioners
care about: if you had predicted function for the proteins that were
unannotated N years ago, how many of your predictions would experiments
have confirmed by now?

## Task

Three horizons. For each, the cohort is the set of proteins (human, yeast,
fly; see paper for construction) that had NO experimentally supported GO
annotation at the horizon date but have at least one today.

| horizon | past GOA releases | cohort |
|---|---|---|
| 2019 | human 192, yeast 93, fly 93 (mid-2019) | 747 |
| 2021 | human 206, yeast 107, fly 107 (mid-2021) | 388 |
| 2023 | human 218, yeast 119, fly 119 (2023-07-12) | 281 |

Submit, for each cohort protein and GO branch, your single most specific
predicted term. Your method may use ANY knowledge that predates the horizon
(sequences, structures, literature, databases as they stood then); using
present-day curated resources is temporal leakage and disqualifies the
leakage-controlled comparison (report such runs separately, as we do for
our domain-augmented arm).

## Files

- cohorts/cohort_<year>.txt - one UniProt accession per line
- truth/truth_<year>.tsv - accession, branch, go_term; the propagated
  current experimental annotation set, frozen at package generation
  (metadata.json records the date and source releases)
- score.py - self-contained scorer (needs go-basic.obo); reports per-branch
  top-1 precision and a permutation p-value
- baselines.md - reference results: a sequence+structure transfer method
  and a past-corpus frequency prior

## Scoring

A prediction hits if the predicted term is in the protein's propagated
truth set (a correct ancestor counts). Because ancestor hits make raw
precision gameable by shallow terms - the frequency prior reaches 0.97 raw
precision on cellular component by predicting a near-root term - the
PRIMARY metric is mean information gain: the information content of the
predicted term (from the horizon's frozen IC table, ic/ic_<year>.tsv) when
correct, zero when wrong, averaged over scored proteins. Shallow correct
guesses earn near-zero bits; wrong specific guesses earn zero. Raw top-1
precision and a permutation p-value are reported alongside. Example:

    python score.py --predictions baselines/ss_2023.tsv \
        --truth truth/truth_2023.tsv --obo go-basic.obo \
        --ic ic/ic_2023.tsv

The baselines/ directory contains our sequence+structure submissions in
the exact expected format.

## Notes and fair-use rules

- Truth is frozen: as GOA accumulates further annotations, scores computed
  against these files remain comparable across submissions. Regenerating
  truth from a newer GOA creates a new benchmark edition; label it as such.
- AlphaFold models postdate some horizons; we treat predicted structure as
  horizon-legal on the grounds that structure is a property of the
  sequence, not knowledge about its function. Disagree? Report with and
  without, as a sensitivity.
- Cohort proteins were unannotated at the horizon partly because they are
  hard to study; absolute precision values are not comparable across
  horizons (cohorts differ), only across methods within a horizon.

Produced from the protein-stack repository; construction details and the
archived-release protocol are described in the accompanying manuscript.
