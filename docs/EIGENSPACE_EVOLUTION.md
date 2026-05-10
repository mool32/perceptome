# Eigenspace evolution — v0.2 (12-PC) → v0.3 (9-PC)

**TL;DR.** The reference eigenspace was rebuilt for perceptome v0.2.0 because the catalog grew from 43 to 44 modules (NPAS4 added). The new eigenspace has 9 PCs satisfying the Kaiser criterion (eigenvalue > 1), versus 12 PCs in v0.1. **PC1 is highly preserved** (`|cos(PC1_v0.1, PC1_v0.3)| = 0.90`), but **PC2-PC3 underwent rotation and reordering** — what was v0.1 PC3 maps best to v0.3 PC2, with mid-range cosines (0.4-0.7). Beyond PC3, the two eigenspaces diverge as they describe different parts of the variance spectrum.

**Practical consequence.** PC1 interpretations port across versions. PC2-PC4 interpretations need re-examination. Numerical magnitudes will differ.

This page exists so that readers comparing the perceptome v0.1 / Paper 3 results to v0.2 tool output (or vice versa) understand the relationship.

---

## What changed

| | v0.1 (Paper 3 release, 2026-05-05) | v0.2 (this release, 2026-05-10) |
|---|---|---|
| Modules | 43 | 44 (+ NPAS4) |
| Cell-type reference | 154 HPA | 154 HPA (same source) |
| Kaiser n_PCs | 12 | 9 |
| Bootstrap stable PCs | PC1-PC3 | PC1-PC6 |
| Top eigenvalue | 9.50 | 12.83 |
| Cumulative variance, top 3 PCs | ~36% | ~51% |
| Spectral fit α (power-law slope, log-log) | (paper3 reports) | 1.46 (R²=0.93) |

## Why fewer PCs?

Two reasons combine:

1. **NPAS4 strengthens an existing axis** rather than introducing a new orthogonal one. NPAS4 readiness in HPA is dominated by the neuronal cell-type class, which already contributes to the same axis as other neuron-relevant modules in v0.1. Adding NPAS4 raises that axis's eigenvalue but does not create a new component.

2. **Slight scoring methodology refinement.** v0.2 uses `mean log1p(nCPM)` per cell type at the catalog-aggregation step (matching v1.4 perceptivity metric design), which is more stable across cell types than the raw per-cell mean used in v0.1's eigenspace construction. This makes intra-module gene panels behave more coherently → tighter correlation block → larger top eigenvalues → fewer PCs need to clear the Kaiser cutoff.

Net: the same biological structure, expressed more compactly. PC1-PC6 of v0.3 are now bootstrap-stable (versus only PC1-PC3 in v0.1), so users get **more reliable** structure across the first six axes despite fewer total PCs.

## Mapping between v0.1 and v0.3 PCs (computed)

For each v0.1 PC, the closest-cosine v0.3 PC. Computed by restricting v0.3 loadings to the 43 v0.1 modules (dropping NPAS4 from v0.3 column), normalizing, and taking argmax|cos| per row. Numbers from `scripts/05_compare_eigenspaces.py`.

| v0.1 PC | Closest v0.3 PC | \|cos\| | Note |
|---|---|---|---|
| PC1 | PC1 | **0.903** | perception-breadth axis well preserved |
| PC2 | PC2 | 0.553 | rotation present |
| PC3 | PC2 | 0.476 | v0.1 PC3 component absorbed into v0.3 PC2 |
| PC4 (cancer convergence axis) | PC4 | **0.708** | cancer convergence axis substantially preserved |
| PC5 | PC5 | 0.297 | weakly mapped |
| PC6 | PC8 | 0.352 | reordered |
| PC7 | PC7 | 0.297 | weakly mapped |
| PC8 | PC9 | 0.528 | reordered |
| PC9 | PC4 | 0.262 | overlap with PC4 axis |
| PC10-PC12 | various | <0.45 | below v0.3 Kaiser cutoff anyway |

**Interpretation.** PC1 (perception breadth) is the most preserved axis — its biological meaning (myeloid+ to germ−) carries directly. PC4 (cancer convergence axis from Paper 3 §3.2) is preserved at 0.71, meaning the cancer-attractor finding ports between versions with modest rotation. PC2-PC3 of v0.1 reorganize into a single dominant axis in v0.3 (likely because the smoother v1.4 scoring procedure causes some modules to co-vary more tightly). PC5+ are not stable enough in either version to expect cross-version mapping.

**For Paper 3 readers:** the perception-breadth ranking on PC1 ports directly. The cancer convergence story on PC4 ports with a 0.71 cosine — qualitative direction preserved, exact loadings differ. PC2/PC3 require re-examination if you want to relate Paper 3 PC2/PC3 narrative to v0.3 output.

## What this means in practice

**For users running the v0.2 tool on new data:** you get the v0.3 9-PC space. PC1 = perception breadth (preserved from Paper 3). PC2-PC3 are dominant stress/identity axes but **rotated** vs Paper 3 numbering. PC4 carries the cancer convergence direction. PC1-PC6 are bootstrap-stable.

**For readers of Paper 3 (v0.1 / 12-PC):** PC1-axis findings reproduce with high fidelity in v0.3. PC4 cancer convergence reproduces with modest rotation. PC2-PC3 narrative needs reinterpretation in v0.3 because of axis reordering. Numerical magnitudes differ across versions in all cases — for direct cross-version comparison, install `perceptome==0.1.0` explicitly.

**For Paper 4.2 (cancer attractor capacity-direction):** the attractor direction is bundled in `pct.reference.attractor` as both per-module Δ (44 modules, framework-stable) and as the v0.3 eigenspace projection (9 PCs). Use `mode='modules'` in `pct.attractor_alignment()` for cross-version consistency, or `mode='eigenspace'` for the v0.3-native projection.

## Backward compatibility

- The v0.1 12-PC `reference_eigenspace.json` is preserved unchanged in `tool/perceptome/perceptome/data/` (the v0.1 release tree). Anyone reproducing Paper 3 numbers should `pip install` v0.1.0 specifically.
- The v0.2 `pct.project()` uses the v0.3 9-PC reference by default; pass `reference=path_to_v01_json` to project against the legacy 12-PC space.
- Disease (RA/AD/IPF/DKD) and aging (inflammaging/collapse) reference vectors from v0.1 are pinned as `*_v02pinned.json` — they are in v0.2 12-PC space and incompatible with v0.3 9-PC. Recompute scheduled for v0.2.1 (raw scRNA-seq data confirmed available locally; rerun planned).

## Why we did NOT keep both eigenspaces as parallel APIs

Having both `pct.project_v01()` and `pct.project_v02()` in the same tool would invite users to mix versions silently — projecting their data into one space, comparing against references in another, and producing meaningless cosines. Single canonical eigenspace per tool version. Users who need v0.1 reproducibility install v0.1 explicitly.

## References

- Paper 3 manuscript (in preparation): perceptome v0.1, 43 modules, 12-PC eigenspace
- Paper 4 (Memory consolidation, in preparation): NPAS4 candidate 44th module, full module composition in Table S4
- v0.3 build script: `tool/perceptome2/scripts/03_build_eigenspace.py`
- v0.3 reference: `tool/perceptome2/perceptome/eigenspace/data/reference_v03.json`
