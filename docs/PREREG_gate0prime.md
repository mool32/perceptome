# Pre-registration: gate 0′ — co-retention of modules beyond the history of their genes

Status: **locked** (`docs/PREREG_LOCK_gate0prime.txt`). Script: `scripts/19_gate0prime_coretention.py`.
Synthetic checks: `tests/test_gate0prime.py`.

## 1. Why gate 0 is replaced

Gate 0 (instruction v2) passed in 3/3 independent curations (C2–C4), but its criterion was nearly
tautological. In Paper 7's P5 the null is almost the same for every module (mean 7.8–8.8, SD 0.90–1.06
in C4), so the module z is essentially the number of species in which the module is intact
(r = 0.98 with `observed` in C4), and the module age is read from the same presence row. "Ancient
modules have higher z" therefore restates "modules found in distant species are found in more
species". Gate 0's PASS shows that this property is curator-independent; it says nothing about modules
being evolutionary units. Recorded as such; not evidence for two layers.

## 2. Question

Are the genes of a module gained and lost **together** more than genes with the same individual
histories (Wagner-type units)? Two-layer claim of Paper 7, framing v3: ancient modules (≥ 1500 Mya)
are coherent units; young modules (≤ 800 Mya) are anti-coherent, assembled from components with
different histories.

## 3. Statistic

- Presence matrix: gene × 15 species, 1 = ortholog by the Paper 7 step03/step04 rule.
- Module coverage per species = fraction of its genes present; module intact if coverage ≥ 0.6;
  module age = largest divergence among species where it is intact (as in the pipeline).
- **S = variance over species of the module's coverage.** With each gene's history fixed, S rises
  when genes co-occur and falls when they replace each other. No threshold is involved, so the sign
  does not depend on whether coverage is above or below K.
- **Null:** each gene replaced by a random gene of the **same gene age** (deepest species with an
  ortholog) **and the same breadth** (number of species with an ortholog); no gene twice in one null
  set; 10,000 sets. z_coh = Φ⁻¹(mid-p quantile of S among null sets).
- Pool: union of all catalog genes **plus ≥ 2,000 background protein-coding genes** run through the
  same pipeline. With smaller pools the exact (age, breadth) classes are too thin and widening them
  biases z_coh downward (shown in synthetic worlds: mean z −0.26 with 760 genes, ≈ 0 with ≥ 2,000).

## 4. Built-in checks per catalog (a catalog failing either is not evaluable)

- widened null classes for > 10 % of its genes;
- calibration: pseudo-modules (each real module's genes replaced by random genes of the same gene age,
  3 per module) must have |mean z_coh| ≤ 0.20.

The pipeline's intact counts and ages, recomputed from the exported presence matrix, must equal
`per_module_Ck.tsv` for every module; otherwise the run stops (`STOP_INPUT_MISMATCH`).

## 5. Criteria (per catalog; ≥ 5 modules in each age group)

- L1 ancient coherent: Wilcoxon signed-rank, z_coh of ancient modules > 0, one-sided p < 0.05.
- L2 young anti-coherent: z_coh of young modules < 0, one-sided p < 0.05.
- L3 ancient > young: Mann–Whitney one-sided p < 0.05.

| verdict | rule over C2, C3, C4 (C0 reported, not counted) |
|---|---|
| **PASS_TWO_LAYERS** | ≥ 2 catalogs with L1, L2 and L3 |
| **PASS_ANCIENT_UNITS_ONLY** | otherwise, ≥ 2 catalogs with L1 and L3 |
| **FAIL** | otherwise |
| **NOT_EVALUABLE** | < 2 evaluable catalogs |

Descriptive: z_coh of HSF1-like modules, Spearman ρ(age, z_coh), the intact-count z with the same
matched null (`z_intact_matched`).

## 6. Validation (tests)

Tautology world (genes of similar origin, independent losses): the old gate-0 criterion passes in
every catalog, gate 0′ never passes and L2 fires in ≤ 2 of 15 catalog tests. Two-layer world (ancient
modules lose genes together, young modules are built from genes lost in complementary species):
PASS_TWO_LAYERS. Ancient-only world: PASS_ANCIENT_UNITS_ONLY.

## 7. Meaning of outcomes

PASS_TWO_LAYERS: Paper 7's v3 framing survives a non-circular test and curator change.
PASS_ANCIENT_UNITS_ONLY: ancient state-monitoring modules are units; "young = borrowed" is unsupported.
FAIL: modules are not evolutionary units beyond the ages of their genes; the evolutionary line is closed.
