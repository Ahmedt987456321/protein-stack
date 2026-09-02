# Baseline results

Primary metric: mean information gain (bits) = IC of the predicted term when correct, 0 when wrong, averaged over scored proteins; IC tables frozen per horizon from the past corpus. Raw top-1 precision shown for interpretability - note the prior's high raw precision but near-zero information gain, which is why raw precision is not the ranking metric.

| horizon | branch | cohort | SS n | SS top-1 | SS gain | prior term | prior top-1 | prior gain |
|---|---|---|---|---|---|---|---|---|
| 2019 | molecular_function | 747 | 305 | 0.351 | 1.208 | GO:0005488 (binding) | 0.682 | 0.557 |
| 2019 | biological_process | 747 | 162 | 0.198 | 1.589 | GO:0009987 (cellular process) | 0.760 | 0.372 |
| 2019 | cellular_component | 747 | 206 | 0.345 | 1.394 | GO:0110165 (cellular anatomical structure) | 0.972 | 0.226 |
| 2021 | molecular_function | 388 | 121 | 0.364 | 2.100 | GO:0005488 (binding) | 0.481 | 0.319 |
| 2021 | biological_process | 388 | 108 | 0.222 | 1.761 | GO:0009987 (cellular process) | 0.804 | 0.553 |
| 2021 | cellular_component | 388 | 129 | 0.372 | 1.582 | GO:0110165 (cellular anatomical structure) | 0.978 | 0.262 |
| 2023 | molecular_function | 281 | 86 | 0.395 | 2.553 | GO:0005488 (binding) | 0.518 | 0.345 |
| 2023 | biological_process | 281 | 77 | 0.221 | 1.785 | GO:0009987 (cellular process) | 0.823 | 0.536 |
| 2023 | cellular_component | 281 | 97 | 0.392 | 1.602 | GO:0110165 (cellular anatomical structure) | 0.977 | 0.268 |
