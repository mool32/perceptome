# Loop layer (v0.4)

**Status: curated, not yet validated.** Nothing in this layer is used by any
scoring function yet; it is the data foundation for planned loop-integrity
metrics. Tiers are expected from cell-line biology and have not been checked
against a Perturb-seq perturbation list.

## Why a separate layer

Every perceptome module is defined by four criteria, one of which is a
negative-feedback element. Until v0.4 that element was recorded
(`feedback_genes`) but never measured, and it could not be measured cleanly:

- only 17 of 44 modules had feedback genes;
- in 15 of those 17, feedback genes were also activity genes, so any
  "targets vs feedback" statistic would correlate a gene with itself;
- several lists mixed generic or misassigned genes (e.g. UBE2I and PIAS4
  under NF-κB, NKD1/2 and BIRC2/3 under Hippo, GADD45A/B/G in two modules,
  `PHD2` as an alias duplicate of EGLN1).

The v0.3 catalog is left untouched so that every published score reproduces;
the loop layer lives in `perceptome/catalog/data/loops_v04.json`, built by
`scripts/12_build_loops_v04.py` (the curation, with references, is in that
script).

## Fields

| field | meaning |
|---|---|
| `feedback_core` | 1–7 negative-feedback genes induced by the module itself |
| `loop_targets` | module output genes, disjoint from `feedback_core` |
| `feedback_mechanism` | how the loop closes |
| `feedback_is_output` | the output *is* the feedback (chaperone titration, clock, ligand catabolism) — target/feedback coupling is undefined, only noise compression applies |
| `feedback_readable` | `false` when feedback is post-translational only (NRF2/KEAP1, mTOR) |
| `perturbations.break` | remove the feedback → output ↑, cell-to-cell noise ↑, coupling lost |
| `perturbations.input_down` | remove the input → output ↓ |
| `perturbations.input_up` | raise the input with the loop intact → output ↑, noise compressed |
| `tier` | A: active in K562/RPE1 with clean perturbations · B: specific context · C: needs stimulus · D: needs ligand/lineage · E: no clean transcriptional loop |

Tier counts: A 8 · B 4 · C 8 · D 8 · E 16.

## API

```python
import perceptome as pct

pct.list_loops("A")            # ['ERK/MAPK', 'HIF', 'HSF1', 'JAK-STAT', 'NRF2', 'SREBP', 'UPR-ATF6', 'UPR-PERK']
pct.get_loop("p53")            # full entry
pct.loop_gene_sets("p53")      # {'feedback': [...], 'targets': [...]}  — disjoint
pct.loop_output_genes("HSF1")  # readout genes; includes feedback when feedback_is_output
pct.loop_overlaps()            # genes shared between modules' targets (biology, not errors)
pct.validate_loops()           # [] when clean
```

## Built-in specificity tests

Some input nodes are shared on purpose and give free specificity checks:
`MBTPS1/2` (UPR-ATF6 and SREBP), `ARNT` (HIF and AhR). A loop-integrity metric
should move both modules when a shared node is knocked down and only one when a
module-specific node is.

## Evidence from Replogle et al. 2022 (metadata only)

`data/evidence/` holds the per-screen testability table and every loop
perturbation in K562_gwps with cell counts and the authors' knockdown
efficiency. Summary: the essential-gene screens (K562_essential, RPE1) test
almost nothing (1 and 2 loops); K562_gwps tests 10 of 28 loops. Adding a
>= 50 % knockdown requirement leaves six loops with a break-vs-input_down
contrast (ERK/MAPK, HIF, JAK-STAT, SREBP, UPR-IRE1, UPR-PERK); p53 is formally
testable in K562_gwps but K562 lacks functional p53, so it serves as a
negative control. The pre-registration built on this is
`docs/PREREG_loop_integrity_v1.md`.

The `tier` field still records the biology-based expectation; it was not
rewritten from the screen.

## Checking a Perturb-seq screen

`scripts/13_check_perturbseq_overlap.py` reads only the metadata of a screen
(e.g. Replogle et al. 2022, `K562_essential` / `K562_gwps` / `rpe1`) and reports,
per loop, which perturbations are present with enough cells and how many
targets and feedback genes are measured.
