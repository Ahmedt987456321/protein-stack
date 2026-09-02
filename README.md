# protein-stack

Using AlphaFold structures to work out what unknown proteins do.

Sequence search stops being reliable once two proteins share less than about
30% of their sequence. Structure search keeps working there, because shape
survives longer than sequence. So this leans on structure exactly where
sequence gives out, mixes in domain and interaction evidence, and keeps every
prediction tagged with how much to trust it and where it came from.

Full numbers, tables, and per-gate detail live under [`results/`](results/) and
[`DECISIONS.md`](DECISIONS.md). The short version:

## Results

Tested on 9,162 human, yeast, and fly proteins with real experimental
annotations (never electronic ones):

- Below 30% identity, adding structure beats sequence alone (+0.086 Fmax).
- The full fusion beats every single kind of evidence (0.61 vs 0.46).
- It produced 27,992 ranked guesses for 2,564 proteins nobody has annotated.
- Time-travel check: built from the 2023 databases, it predicts annotations
  that only appeared by 2026 at 0.40 top-1 precision (p = 0.001).

Rerunning the whole thing on six species across four kingdoms (42,000 proteins)
holds every one of these.

Five smaller studies reuse the same structures and graph: a pseudo-enzyme
catalog, a measurement of how much curated databases inflate retrospective
tests, pathogenic variants sitting at predicted protein interfaces (23% vs 13%
for benign controls, p = 0.01), pairs that interact both genetically and
physically, and protein flexibility teased out of AlphaFold. Each is written up
under [`results/`](results/). None is confirmed in a lab.

## Running it

```bash
pip install -r requirements.txt
python run_all.py
```

Nineteen steps, one command, one machine; downloads are cached so a rerun picks
up where it stopped. It needs MMseqs2, Foldseek, and fpocket (plus localColabFold
for the complex folding) - see [`docs/WSL_GPU_NOTES.md`](docs/WSL_GPU_NOTES.md)
for the Windows/WSL setup. Settings are in `config.yaml`. Run the tests with
`python tests/test_core.py`.

## One thing to be clear about

Everything here is a prediction, labelled as one wherever it appears. Ground
truth is experimental evidence only, and nothing has been checked in a lab -
each result is a lead for someone to test, not a finding. Falsified guesses stay
in the graph next to the rule that made them.

Inputs are all public: UniProt, AlphaFold DB, GOA, the Gene Ontology, InterPro,
STRING, BioGRID, ClinVar. License: MIT ([`CITATION.cff`](CITATION.cff)).
