# Structural template validation of interactions

For each STRING pair whose partners both have an AlphaFold model, we test whether the two models match different chains of the same experimental PDB complex (TM >= 0.5). This is a GPU-free proxy for predicting the dimer: a shared-complex template both corroborates the interaction and names an interface.

| set | pairs | with structural template | fraction |
|---|---|---|---|
| STRING (high-confidence) | 480311 | 144472 | 30.1% |
| random control (same proteins) | 20000 | 82 | 0.4% |

**Enrichment of real interactions over random: 73.4x.** If large, structure independently supports the STRING network and supplies candidate interfaces.

Interactions involving a dark (unannotated) protein that have structural template support: **8884** - novel, structurally backed interaction hypotheses. Top 15 by confidence:

| protein A | protein B | template PDB | min TM |
|---|---|---|---|
| P57784 | Q62189 | 3jb9 | 1.00 |
| P61965 | Q9D7H2 | 2co0 | 1.00 |
| P01594 | P0DP04 | 6axl | 1.00 |
| P00417 | Q9VAW7 | 8ugj | 1.00 |
| P38060 | Q8JZS7 | 3mp5 | 1.00 |
| Q62189 | Q9CQI7 | 7dwh | 1.00 |
| D3YZG8 | P18155 | 6jid | 1.00 |
| P35561 | P52187 | 1u4f | 1.00 |
| P52187 | P52189 | 1u4f | 1.00 |
| O22795 | P56806 | 5mmm | 1.00 |
| P56804 | P56806 | 5mmm | 1.00 |
| P62305 | P62315 | 4pjo | 1.00 |
| P62305 | P62320 | 4pjo | 1.00 |
| P18935 | Q9VAW7 | 8ugj | 1.00 |
| Q9VAW7 | Q9VXI6 | 8ugj | 1.00 |

Caveat: template-based, not de novo AlphaFold-Multimer; a shared fold-pair template is evidence of feasibility and a likely interface, not proof of a specific cellular complex.

## De-novo prediction on the GPU (AlphaFold2-Multimer)

The 15 smallest template-supported dark-protein interactions were folded as
complexes with AlphaFold2-Multimer (ColabFold, MMseqs2-API alignments) on a
laptop RTX 3050 Ti (4 GB). All 15 completed; predictions ranked by interface
confidence (ipTM):

| pair | dark partner | combined len | template PDB | ipTM | pTM | mean pLDDT |
|---|---|---|---|---|---|---|
| P62309-P62322 | P62309 | 167 | 6v4x | 0.85 | 0.86 | 90.1 |
| P62305-P62309 | P62305 | 168 | 4pjo | 0.84 | 0.85 | 90.0 |
| P62307-P62309 | P62307 | 162 | 8i0t | 0.69 | 0.77 | 88.1 |
| Q9V998-Q9VJ33 | Q9V998 | 157 | 6njd | 0.43 | 0.62 | 78.4 |
| Q39236-Q9FLM8 | Q39236 | 157 | 7lbm | 0.42 | 0.58 | 82.1 |
| P0A978-P0A982 | P0A982 | 140 | 3i2z | 0.37 | 0.60 | 79.7 |
| (9 more) | | | | 0.07-0.24 | | |

**3 of 15 reached ipTM >= 0.5** (confident interface). Honest reading:

- The three confident predictions (P62305/P62307/P62309/P62322) are Sm-class
  small nuclear ribonucleoprotein subunits (SNRPE/F/G family), which are
  well characterised in the literature as forming the heptameric Sm ring.
  They are "dark" here only in the narrow sense that our sampled GOA snapshot
  lacked EXPERIMENTAL evidence codes for them - they are not biologically
  unknown. AF2's high confidence is therefore method validation (the pipeline
  confidently recovers real complexes) rather than novel discovery.
- The genuinely obscure candidates (Arabidopsis Q9SYA6, Drosophila Q9VC49,
  etc.) mostly returned low ipTM, i.e. AF2 would not confidently predict a
  specific complex for them despite STRING and a PDB fold-template.
- Secondary honest finding: a PDB fold-template (both partners match chains
  of one complex) is necessary-ish but NOT sufficient for AF2 interface
  confidence - 12 of 15 template-supported pairs still scored ipTM < 0.5.
  Template support and de-novo interface confidence are different, weaker and
  stronger, forms of evidence.

Scope: 15 small complexes on a 4 GB card, template-supported subset only -
a proof-of-concept, not a proteome-wide screen. The pipeline (structural
screen -> shortlist -> de-novo AF2-Multimer) is what generalises.
