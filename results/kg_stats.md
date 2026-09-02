# Knowledge graph - data/kg.db

2467203 edges, 122412 nodes.

| predicate | evidence | edges |
|---|---|---:|
| structure_similar_to | COMPUTATIONAL | 750928 |
| interacts_with | COMPUTATIONAL | 525761 |
| has_function | EXPERIMENTAL | 435942 |
| sequence_similar_to | COMPUTATIONAL | 365730 |
| has_domain | COMPUTATIONAL | 216441 |
| predicted_function | HYPOTHESIS | 69360 |
| electronic_function | COMPUTATIONAL | 39262 |
| domain_implies_function | CURATED | 30122 |
| hypothesis_unconfirmed | COMPUTATIONAL | 23063 |
| has_druggable_pocket | COMPUTATIONAL | 7046 |
| hypothesis_supported | EXPERIMENTAL | 3538 |
| hypothesis_contradicted | EXPERIMENTAL | 10 |

## Demo: top structure-only hypothesis

`A6NML5` predicted_function `GO:0005515` (score 0.50) - supported only by structural neighbours:

- A6NML5 **electronic_function** GO:0016020  (COMPUTATIONAL, GOA-IEA, 1.00)
- A6NML5 **has_druggable_pocket** pocket  (COMPUTATIONAL, fpocket, 0.75)
- A6NML5 **interacts_with** O95876  (COMPUTATIONAL, STRING, 0.78)
- A6NML5 **predicted_function** GO:0005515  (HYPOTHESIS, fusion-v0.2, 0.50)
- A6NML5 **predicted_function** GO:0005886  (HYPOTHESIS, fusion-v0.2, 0.50)
- A6NML5 **predicted_function** GO:0035579  (HYPOTHESIS, fusion-v0.2, 0.50)
- A6NML5 **predicted_function** GO:0005802  (HYPOTHESIS, fusion-v0.2, 0.43)
- A6NML5 **structure_similar_to** Q96HJ5  (COMPUTATIONAL, Foldseek, 0.71)
- A6NML5 **structure_similar_to** Q9H2W1  (COMPUTATIONAL, Foldseek, 0.61)
- A6NML5 **structure_similar_to** Q969K7  (COMPUTATIONAL, Foldseek, 0.59)
