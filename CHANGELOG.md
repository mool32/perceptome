# Changelog

All notable changes to perceptome will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [0.2.0] — 2026-05-10

Major restructure. Adds the perceptivity / capacity layer (factor 1 of the two-factor framework), validity scorecard, attractor capacity-direction reference, and a focused 9-anchor drug operation. Replaces the v0.1 flat module with per-concern subpackages.

### Added

- **44-module catalog** (was 43). NPAS4 added as 44th canonical module per Paper 4 Table S4 (neuron-perception, sensor Ca²⁺/CaMKII, TF NPAS4, 14 activity targets).
- **`pct.perceptivity` subpackage** — 5-vector capacity / engagement / specialization metric per cell type × module:
  - `R_continuous`, `A_continuous`, `C_continuous = R − A`, `headroom = A_max(M) − A`
  - `BS` (within-class breadth), `GS` (global breadth), `I` (mean intensity per engaged module), `spec_quadrant`
  - `capacity_floor()` classifier (saturated_blocked_up | capacious | intermediate | no_data) — upward-asymmetric per Paper 4.8 holdout PARTIAL.
  - `predict_engagement()` — two-factor framework, factor 1 (capacity) only; factor 2 (operation intensity) reserved for context-specific modeling outside the tool.
- **HPA perceptivity reference** (154 × 44 R/A/C/headroom matrices), precomputed and bundled.
- **`pct.validity` subpackage** — random-200 / housekeeping / cell-cycle null scorecard for perturbation-vs-control comparisons. Standard practice since Paper 4.5 v1.2 amendment caught a baseline-shift artifact in adult-intestinal-epithelium analysis.
- **`pct.reference.attractor`** — Cancer capacity-direction (Paper 4.2 P3 PASS, Sun 2021 paired HCC cohort 4/6 cell-types). Locked attractor cluster (8 cells), per-module Δ vector, eigenspace projection.
- **`pct.drugs` subpackage** — 9-anchor reference table + `activity_layer_screen()` operation. The 9 anchors are: MEKi/ERK, Proteasomei/UPR-PERK, CDKi/Cell Cycle, EGFRi/ERK, PI3Ki/PI3K-PTEN, IKKi/NF-κB (6 validated rescues) + HSP90i/HSF1, HSP90i/UPR-ATF6, HSP90i/NRF2 (3 stress-axis controls). Procedure mirrors Block 5 v1.2 of Paper 4.1 (1000 random non-panel drugs as background, BH FDR q<0.10).
- **9-PC v0.3 eigenspace** rebuilt from 154 × 44 R matrix (NPAS4 included). PC1-PC6 stable per 100 bootstraps, PC7+ probable. Top 5 eigenvalues: 12.83, 6.04, 3.42, 2.51, 2.14. Spectral fit α=1.46, R²=0.93.
- **42 unit tests** covering catalog, scoring, perceptivity (with paper4.4/4.5 a-priori predictions reproduced), eigenspace, attractor, validity, drugs.

### Changed

- **Restructured** from single flat `perceptome` module into per-concern subpackages: `catalog`, `score`, `perceptivity`, `eigenspace`, `compare`, `reference`, `validity`, `drugs`, `network`. Top-level namespace remains flat (`pct.score_modules`, `pct.predict_engagement`, `pct.activity_layer_screen` etc.) — no breaking change for users who import via the flat API.
- `compare_to_references()` no longer accepts `include_drugs=True` — the cosine drug matching it would have done was falsified by Paper 4.1. Use `pct.activity_layer_screen()` and `pct.drug_anchors()` for the validated drug operation.
- `infrastructure_regime()` now includes NPAS4 in the perception module list.

### Removed

- **19,804 CMap drug vectors** (`drug_vectors.csv`). The pre-computed projection of all CMap compounds into the eigenspace was a pre-computed asset for drug-disease cosine matching, an operation falsified by Paper 4.1 in three independent formulations (cosine v1.0, subspace v1.1, module overlap v1.2). Replaced by the 9-anchor `drug_anchors` table + activity-layer screen operation. The original asset is preserved in the v0.1.0 source tree (`tool/perceptome/perceptome/data/drug_vectors.csv`) for reproducibility of paper4.1 evidence trail.

### Recomputed reference vectors (NEW in v0.2.0)

- **Disease vectors (RA/AD/IPF/DKD)** rebuilt on v0.3 eigenspace from raw scRNA-seq:
  - RA: CELLxGENE PBMC (Zhang et al.), 108K cells, 6 cell-type groups
  - AD: Grubman 2019 GSE138852 entorhinal, 13K cells, 5 cell-type groups
  - IPF: Habermann 2020, 89K cells, 9 cell-type groups
  - DKD: CELLxGENE diabetic kidney, 39K cells, 9 cell-type groups
- **Aging reference** rebuilt from CELLxGENE Census full_blood + full_bone_marrow with `age_bin` annotation. Inflammaging direction (blood old−young) and collapse direction (bone marrow old−young) projected into v0.3 eigenspace.
- **Biological refinement vs v0.1**: shared aging modules across blood + bone marrow narrowed from 4 (v0.1: UPR-ATF6, ERK/MAPK, NF-κB, GR) to 2 (v0.3: ERK/MAPK, GR). Reason — UPR-ATF6 actually goes opposite directions in the two tissues (blood +0.054 inflammatory, bone marrow −0.279 collapse), which the v0.3 mean_raw scoring on raw data caught. ERK/MAPK + GR are the genuinely shared perception-axis aging modules; UPR-ATF6 is tissue-specific.
- Build scripts: `scripts/06_recompute_disease_vectors.py`, `scripts/07_recompute_aging_reference.py`. Reproducible from raw h5ad data on disk.

### Honest scope notes (for documentation discipline)

- Drug perturbation analysis: narrow validated scope (9 anchors, activity-layer screen). Do **not** use the tool for drug-disease cosine matching — see Paper 4.1 §2.1-2.3 for three independent falsifications. Five additional theoretical framings (timescale, TF-autoregulation, equilibrium framing, FDA prediction, dynamics blindness) were tested on 2026-05-09 and all KILLED. See `PAPER4_1_CANONICAL_RESULTS.md` for the full audit trail (`killed_hypotheses.csv` for the per-row record).
- Capacity-floor predictor: scope is upward saturation only. Active downward suppression by specific signaling perturbations (e.g. atRA → UPR-ATF6) IS allowed and observed (Paper 4.8 holdout PARTIAL); the predictor does not say cells with saturated baselines cannot move at all, only that they cannot ramp UP further.
- Substrate-series scope: architecture engagement observed in active differentiation / acute remodeling contexts (neurons, plasmablast, Th1, muscle hypertrophy, organoid regeneration). NOT observed in steady-state homeostatic terminally-differentiated cells (adult intestinal goblet, Paneth, enterocyte; Paper 4.5 closed FAIL 2026-05-10).

## [0.1.0] — 2026-05-05

First public release. Accompanies Paper 3 ("The perceptual eigenspace of the cell: 43 transcriptional modules organize cancer and aging").

### Added

- 43-module catalog with core, activity, and feedback gene lists per module
- Three scoring methods: `mean_raw` (recommended), `mean_zscore`, `scanpy` (score_genes wrapper)
- Canonical 12-PC eigenspace derived from 154 Human Protein Atlas cell types
- Eigenspace projection for arbitrary RNA-seq input via `pct.project()`
- Pre-computed reference vectors:
  - Cancer (lung, HCC, BRCA matched cell-type pairs)
  - Disease (AD, DKD, RA, IPF perturbation vectors)
  - Aging (blood, bone marrow tissue-specific aging vectors)
  - Drug pharmacology (19,804 CMap/LINCS L1000 compound projections) — *removed in v0.2 per Paper 4.1 falsification*
- Module co-variation network analysis with rich-club detection
- Infrastructure regime classification: 4-regime (supply chain / firefighting / collapse / unsupported)
- Distribution-shape metrics per cell type: Shannon entropy, Gini coefficient, kurtosis, skewness
- Module Importance Index (MII) and Module Discriminative Value (MDV) computation
- Plotting helpers: eigenspace scatter, trajectory, heatmap, regime distribution
- 26 unit tests covering scoring determinism, gene-name mapping, eigenspace consistency, reference loading
