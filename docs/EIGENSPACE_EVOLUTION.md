# Eigenspace evolution — first Paper 3 draft (43 modules, 12 PCs) → final Paper 3 / v0.3 reference (44 modules, 9 PCs)

**TL;DR.** The first draft of Paper 3 ("The perceptual eigenspace of the cell: 43 transcriptional modules organize cancer and aging"; tool v0.1.0) used 43 modules and a 12-PC eigenspace. The final Paper 3 — Spiro (2026), *Diverse cancers converge on one region of normal cell-state space in a cancer-blind coordinate system* ([Zenodo 10.5281/zenodo.20542130](https://doi.org/10.5281/zenodo.20542130)) — and this tool's reference use 44 modules (NPAS4 added) and **9 PCs** satisfying the Kaiser criterion (eigenvalue > 1). **PC1 is highly preserved** between the two (`|cos(PC1_v0.1, PC1_v0.3)| = 0.90`), but **PC2–PC3 underwent rotation and reordering** — what was v0.1 PC3 maps best to v0.3 PC2, with mid-range cosines (0.4–0.7). Beyond PC3 the two eigenspaces diverge.

**Practical consequence.** PC1 interpretations port across versions. PC2–PC4 interpretations need re-examination. Numerical magnitudes differ. And one framing of the first draft does not survive: **the final paper does not describe convergence as an axis** (see "How the final paper describes the axes").

This page exists so that readers comparing the first-draft results to current tool output (or vice versa) understand the relationship.

---

## What changed

| | first draft (v0.1, 2026-05-05) | final Paper 3 / this tool's reference |
|---|---|---|
| Modules | 43 | 44 (+ NPAS4) |
| Cell-type reference | 154 HPA | 154 HPA (same source) |
| Kaiser n_PCs | 12 | 9 |
| Bootstrap stable PCs | PC1–PC3 | PC1–PC6 (0.63–0.99); PC7–PC9 probable (0.40–0.54) |
| Top eigenvalue | 9.50 | 12.83 |
| Cumulative variance, top 3 PCs | ~36% | ~51% |
| Cumulative variance, the nine retained PCs | — | 73.8% |
| Spectral fit α (power-law slope, log-log) | (draft reports) | 1.46 (R²=0.93) |

## Why fewer PCs? What NPAS4 does and does not explain

An earlier version of this page attributed the drop from 12 to 9 components partly to NPAS4 "strengthening an existing axis". The final paper's leave-one-out (Note S2) shows that is not the explanation: **rebuilding the current pipeline on 43 modules (without NPAS4) reproduces all nine components, with a per-PC eigenvector cosine ≥ 0.95 (mean 0.99 over PC1–PC6)** — reproduced by `tests/test_paper_reference.py` (minimum 0.951, mean 0.993). NPAS4 is the single neuron-restricted module; its score is near zero in non-neural cell types and it does not shape the cross-tissue space. All cross-tissue statistics therefore transfer unchanged between the 43- and the 44-module catalog.

What remains is the scoring change: the current reference uses `mean log1p(nCPM)` per cell type at the catalog-aggregation step (matching the v1.4 perceptivity metric design), which is more stable across cell types than the raw per-cell mean used for the first draft's eigenspace construction. That makes intra-module gene panels behave more coherently — a tighter correlation block, larger top eigenvalues, fewer PCs clearing the Kaiser cutoff. The same biological structure, expressed more compactly, with PC1–PC6 now bootstrap-stable rather than PC1–PC3.

## Mapping between v0.1 and v0.3 PCs (computed)

For each v0.1 PC, the closest-cosine v0.3 PC. Computed by restricting v0.3 loadings to the 43 v0.1 modules (dropping NPAS4 from v0.3 column), normalizing, and taking argmax|cos| per row. Numbers from `scripts/05_compare_eigenspaces.py`.

| v0.1 PC | Closest v0.3 PC | \|cos\| | Note |
|---|---|---|---|
| PC1 | PC1 | **0.903** | the leading axis is well preserved |
| PC2 | PC2 | 0.553 | rotation present |
| PC3 | PC2 | 0.476 | v0.1 PC3 component absorbed into v0.3 PC2 |
| PC4 (the draft's "cancer convergence axis") | PC4 | **0.708** | the axis itself is substantially preserved; the draft's reading of it as a convergence axis is not carried into the final paper |
| PC5 | PC5 | 0.297 | weakly mapped |
| PC6 | PC8 | 0.352 | reordered |
| PC7 | PC7 | 0.297 | weakly mapped |
| PC8 | PC9 | 0.528 | reordered |
| PC9 | PC4 | 0.262 | overlap with PC4 axis |
| PC10-PC12 | various | <0.45 | below v0.3 Kaiser cutoff anyway |

PC5+ are not stable enough in either version to expect cross-version mapping.

## How the final paper describes the axes

In the standardized coordinates of the final paper (and of `pct.project`):

| Axis | Share of variance* | Runs from … to … | Independent check |
|---|---|---|---|
| PC1 | 31% | inflammatory cells → terminal germline cells | an external immune-signalling signature orders the axis monotonically (Spearman ρ = −0.89) |
| PC2 | 21% | quiescent → proliferative cells | translational output (ρ = +0.75) |
| PC3 | — | signaling-dominated immune cells → metabolic epithelia | ρ = +0.50 |
| PC5 | — | — | cholesterol signature (ρ = +0.54) |

\*As printed in Fig. 1B. The shipped reference gives 29.2% and 13.7%; the nine-component total (73.8%) matches. The printed labels presumably come from the single-cell-scaled matrix (see [`PAPER_MAP.md`](PAPER_MAP.md)). A biological-family label not used to build the space explains 31–70% of cell position on every one of PC1–PC6 (η², each against its own label-permutation null).

The bundled `pct.project` reproduces the orientation: neutrophils, monocytes and macrophages sit at the negative end of PC1; oocytes and megakaryocytes at the positive end.

**Convergence is not an axis.** The paper states that convergence "is a high-dimensional nearest-neighbor property, not a 2-D location": a raw PC1–PC2 projection recovers only 8.2% core assignment for B-ALL against 60.8% in the full nine dimensions (6.8% against 23.9% pooled over all 25 cancers). The first draft's labels for PC1 ("perception breadth") and PC4 ("cancer convergence axis") are therefore historical; the README figure that plotted PC1 × PC4 has been replaced by the paper's objective plane (cosines to the two core cells, `scripts/09_build_eigenspace_figure.py`).

## What this means in practice

**For users running the tool on new data:** you get the 9-PC space, placed with the standardized projection. PC1–PC6 are bootstrap-stable. Interpret PC1–PC3 as in the table above, and use `pct.convergence` — not a PC — for the cancer question.

**For readers of the first draft (v0.1 / 12-PC):** PC1-axis findings reproduce with high fidelity. PC2–PC3 narrative needs reinterpretation because of axis reordering. Numerical magnitudes differ across versions in all cases — for direct comparison install `perceptome==0.1.0` explicitly.

**For Paper 4.2 (legacy cancer capacity-direction):** the direction is bundled in `pct.reference.attractor` as both a per-module Δ (44 modules, framework-stable) and its eigenspace projection (9 PCs, standardized coordinates since this release). Use `mode='modules'` in `pct.attractor_alignment()` for cross-version consistency, or `mode='eigenspace'` for the native projection.

## Backward compatibility

- The v0.1 12-PC `reference_eigenspace.json` is preserved unchanged in the v0.1 release tree. Anyone reproducing first-draft numbers should `pip install` v0.1.0 specifically.
- `pct.project()` uses the 9-PC reference by default; pass `reference=path_to_v01_json` to project against the legacy 12-PC space.
- **`pct.project()` now standardizes** (module scores z-scored against the reference calibration before the loadings), as in the paper's Methods. Versions ≤ 0.2.2 projected raw scores, which left every cell sharing one large offset (mean pairwise cosine between the 154 normal cell types 0.97, instead of ~0). `standardize=False` restores the old output. The disease, aging and attractor vectors stored with the package were re-expressed in the standardized coordinates from their stored per-module deltas (`scripts/12_restandardize_reference_vectors.py`); cosines with them from older versions are not comparable.
- The disease (RA/AD/IPF/DKD) and aging (inflammaging/collapse) vectors were recomputed from raw scRNA-seq on the 9-PC reference in v0.2.0; the 12-PC v0.1 versions are not shipped.

## Why we did NOT keep both eigenspaces as parallel APIs

Having both `pct.project_v01()` and `pct.project_v02()` in the same tool would invite users to mix versions silently — projecting their data into one space, comparing against references in another, and producing meaningless cosines. Single canonical eigenspace per tool version. Users who need v0.1 reproducibility install v0.1 explicitly.

## References

- Spiro, T. (2026). *Diverse cancers converge on one region of normal cell-state space in a cancer-blind coordinate system.* Zenodo, [10.5281/zenodo.20542130](https://doi.org/10.5281/zenodo.20542130) — final Paper 3: 44 modules, nine-dimensional reference space, Notes S1 (eigenspace), S2 (NPAS4 leave-one-out), Fig. S1 (stability), Fig. S11 / S14 (axes).
- First Paper 3 draft (v0.1.0): 43 modules, 12-PC eigenspace.
- Paper 4 (Memory consolidation, in preparation): NPAS4 candidate 44th module, full module composition in Table S4.
- Build script: `scripts/03_build_eigenspace.py`; reference: `perceptome/eigenspace/data/reference_v03.json`.
