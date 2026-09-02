# v0.3 - dark-proteome sweep - report

Dark set: **7128 proteins** (GOA proteins with zero experimental GO evidence), knowledge base: 42232 annotated proteins, fusion params from v0.2: `{'tau': 0.35, 'alpha': 0.7, 'beta': 0.0, 'gamma': 1.0, 'tau_d': 0.8, 'delta': 0.0, 'val_macro_fmax': 0.6006}`.

## Hypothesis feed

- 69360 hypotheses across 7128 proteins (HIGH 27400, MEDIUM 26522, LOW 15438)
- **Structure-only proteins: 13** - no sequence neighbours, no InterPro domains; AlphaFold+Foldseek is the only evidence. These are the highest-novelty targets.

## Corroboration vs held-aside electronic annotations (top MF hypothesis, n=4544 checkable)

| tier | agree | total | rate |
|---|---:|---:|---:|
| HIGH | 3819 | 4182 | 91.3% |
| MEDIUM | 208 | 299 | 69.6% |
| LOW | 15 | 63 | 23.8% |
| shuffled baseline | - | - | 10.5% |

Electronic annotations are themselves computational; agreement is corroboration, not validation. Disagreements in the HIGH tier are the interesting review queue.

**Gate** (HIGH-tier corroboration >= 2x shuffled baseline, and >= 100 HIGH hypotheses): HIGH 91.3% vs baseline 10.5% (ratio 8.7x), HIGH count 27400 -> **PASSED**
