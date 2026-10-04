# Pre-registration v2: homeostatic locality of knockdown responses

Status: **locked** (see `docs/PREREG_LOCK_v2.txt`). Analysis code: `scripts/17_xatlas_pseudobulk.py`
(pseudobulk builder) and `scripts/18_homeostatic_locality_prereg_v2.py` (test). Synthetic checks:
`tests/test_homeostatic_locality_prereg.py`.

## 1. Where the hypothesis comes from (honest provenance)

Exploration on Replogle 2022 K562_gwps (`docs/EXPLORE_dose_buffering.md`), data already used by the
locked loop-integrity test v1 (NULL). The response of each knockdown was split into the module's own
activity genes versus the whole transcriptome ("locality"). Three groupings of modules were tried on
those same data:

| grouping | adjusted median, group vs rest | one-sided MW p |
|---|---|---|
| curated tier A loops | +0.46 vs −0.29 | 0.003 (circular: tier A was curated with K562 activity in mind) |
| non-empty transcriptional `feedback_core` | +0.24 vs −0.42 | 0.015; **0.070** once knockdowns of the module's own targets are removed |
| catalog category `B_interoceptive` | +0.20 vs −0.42 | 0.026; **0.007** once own-target knockdowns are removed |

The third grouping is the one locked here. It was **chosen after seeing K562**, but the category itself
predates all Perturb-seq work (catalog v0.2/v0.3, defined from biology: modules that sense an internal
state the cell regulates — sterols, oxygen, redox, unfolded protein, heat, energy, nutrients, DNA
damage, cytosolic DNA, Ca²⁺). It also coincides largely with the "ancient state-monitoring layer" of
Paper 7. K562 results are therefore not evidence; only the data below are.

## 2. Hypothesis

**H (homeostatic closure):** a knockdown of the machinery of an interoceptive module changes that
module's own target genes more, relative to expression-matched genes elsewhere in the transcriptome,
than a knockdown of the machinery of a non-interoceptive module (external-signal relays and nuclear
receptors), after adjusting for the baseline expression of the module's targets.

Mechanistic reading: a loop that senses its own output absorbs a perturbation inside the loop, so the
response stays confined to the loop; a relay passes the perturbation on.

## 3. Data

X-Atlas/Orion (Xaira Therapeutics, 2025): genome-wide CRISPRi Perturb-seq, FiCS platform, all
protein-coding genes, two lines — **HCT116 (primary)** and **HEK293T (replication)**. Different lab,
platform and cell lines from Replogle 2022. Files: the per-line filtered dual-guide h5ad (raw UMI),
read lazily (HF / Figshare+). At lock time **no expression value of X-Atlas has been read by anyone in
this project**; the builder's `--inspect` mode, which reads only column names and category labels, is
the only access permitted before the lock.

## 4. Pseudobulk (script 17)

- Knockdowns read: union over the 44 modules of core/sensor/tf/feedback genes plus loop `feedback_core`
  (466 symbols; aliases in the script). No other perturbation is read.
- Controls: seeded sample of 10,000 non-targeting cells (all if fewer). log-CP10k; per batch μ and σ
  (σ floored at 0.05) if the batch has ≥ 200 sampled controls, otherwise pooled.
- Per knockdown: mean over its cells of the per-cell z-score; fold = mean target CP10k in knockdown
  cells / mean in controls; control mean UMI of the target; number of cells.
- Column flags (`--target_col`, `--control_label`, `--batch_col`, `--gene_col`, `--filter_col`) are
  set from `--inspect` output only and recorded in `pseudobulk_meta.json`.

## 5. Unit-level metric (script 18)

- Readout universe: genes with control mean ≥ 0.01 UMI/cell.
- Module targets A: catalog `activity` genes in the universe, minus the knocked-down gene; ≥ 3 needed.
- Efficient knockdown: ≥ 50 cells, target control mean ≥ 0.10 UMI/cell, fold ≤ 0.50.
- **Knockdowns of the module's own activity genes are excluded** (they are the readout, not machinery).
- Per (knockdown, module): L = log2( mean z² over A / median of 500 null means ), each null set drawing,
  for every gene of A, a random gene of the same control-expression decile outside A (seeded by
  hash of module, knockdown). This removes the detectability confound.
- Module score: median L over its efficient knockdowns (≥ 1). Baseline: log10 mean control UMI of A.

## 6. Primary test and decision

Per line: OLS of module score on baseline; one-sided Mann–Whitney of residuals, `B_interoceptive`
greater than all other modules. A line is evaluable with ≥ 8 modules in each group.

| outcome | rule |
|---|---|
| **CONFIRMED** | HCT116 evaluable and p < 0.05, and HEK293T evaluable and p < 0.05 |
| **PARTIAL** | HCT116 p < 0.05; HEK293T not significant or not evaluable |
| **NULL** | HCT116 evaluable, p ≥ 0.05 |
| **NOT_EVALUABLE** | HCT116 not evaluable |

**Negative control N1 (off-module readout):** each knockdown is scored against the activity genes of a
random module it does not belong to and whose activity set does not overlap any activity set of its own
modules. If the same contrast is p < 0.05 in either line, the verdict gets the suffix `_NONSPECIFIC`
(the effect would be about how strongly interoceptive knockdowns move the transcriptome, not about
locality) and cannot be reported as CONFIRMED.

## 7. Secondary analyses (reported, no decision)

S1 `feedback_core` grouping; S2 without baseline adjustment; S3 including own-target knockdowns
(the K562 definition); S4 modules with ≥ 2 knockdowns; S5 fold ≤ 0.60 (HEK293T knockdown is weaker,
median 51.5 %); descriptive Paper 7 cohorts from the per-module table.

## 8. Expected power and validation

- Synthetic worlds (tests): null world never CONFIRMED (0/30; primary p < 0.05 in 3 %), detectability
  confound world 0/10 primary hits, planted effect CONFIRMED.
- Realistic power, bootstrapping K562 module residuals: per line 0.83 (16 vs 20 modules) at the K562
  effect size, 0.54 at half of it; both lines jointly ≈ 0.69 and ≈ 0.3. K562 effect sizes are likely
  inflated (definition chosen on those data), so a NULL here is a weak negative, a CONFIRMED a strong
  positive.

## 9. What changes after the result

CONFIRMED: homeostatic closure becomes the lead finding — a property of signalling architecture visible
only in the shape of the response, not in mean activity; next steps are dose (knockdown strength vs
locality) and the Paper 7 age link. PARTIAL: report as line-dependent. NULL: the K562 pattern is treated
as a product of exploration, and the locality direction is closed.
