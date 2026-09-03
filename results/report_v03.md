# v0.3 - dark-proteome sweep - report

Dark set: **7214 proteins** (GOA proteins with zero experimental GO evidence), knowledge base: 41157 annotated proteins, fusion params from v0.2: `{'tau': 0.35, 'alpha': 0.7, 'beta': 0.0, 'gamma': 1.0, 'tau_d': 0.8, 'delta': 0.0, 'val_macro_fmax': 0.5955}`.

## Hypothesis feed

- 70498 hypotheses across 7214 proteins (HIGH 25500, MEDIUM 27908, LOW 17090)
- **Structure-only proteins: 9** - no sequence neighbours, no InterPro domains; AlphaFold+Foldseek is the only evidence. These are the highest-novelty targets.

## Corroboration vs held-aside electronic annotations (top MF hypothesis, n=4679 checkable)

| tier | agree | total | rate |
|---|---:|---:|---:|
| HIGH | 4002 | 4378 | 91.4% |
| MEDIUM | 153 | 242 | 63.2% |
| LOW | 9 | 59 | 15.3% |
| shuffled baseline | - | - | 11.1% |

Electronic annotations are themselves computational; agreement is corroboration, not validation. Disagreements in the HIGH tier are the interesting review queue.

**Gate** (HIGH-tier corroboration >= 2x shuffled baseline, and >= 100 HIGH hypotheses): HIGH 91.4% vs baseline 11.1% (ratio 8.2x), HIGH count 25500 -> **PASSED**
