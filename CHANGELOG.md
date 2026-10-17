# Changelog

All notable changes to perceptome will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/), and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

Alignment with the final Paper 3 — Spiro (2026), *Diverse cancers converge on one region of normal cell-state space in a cancer-blind coordinate system* ([Zenodo 10.5281/zenodo.20542130](https://doi.org/10.5281/zenodo.20542130), manuscript v11) — and correction of statements it retracts. [`docs/PAPER_MAP.md`](docs/PAPER_MAP.md) maps every part of the paper onto the code and lists which printed numbers the bundled data do and do not reproduce.

### Changed — behaviour (read before upgrading)
- **`pct.project` standardizes by default.** The paper places a cell in the space by z-scoring each module against the reference calibration and projecting through the fixed loadings. Versions ≤ 0.2.2 multiplied raw scores by the loadings, leaving all cells with one large shared offset (mean pairwise cosine between the 154 normal cell types 0.97; 0.002 once standardized, matching the paper's near-orthogonal states), labelling every HPA cell type `low_perception`, and reversing the paper's PC1 orientation. `project(..., standardize=False)` restores the old output. Modules missing from the input are now neutral (placed at the reference mean) rather than raw zero; `missing_modules=` marks uncovered modules.
- `project` warns when per-cell-scale scores meet the pseudobulk calibration (paper Notes S-B, S-H); `n_pcs` defaults to all PCs of the reference.
- **`compare_conditions()["delta_eigenspace"]`** is computed with `pct.project_shift` (module delta divided by the reference SD, then the loadings), so it equals the difference between the projected condition means.
- **The bundled attractor / disease / aging eigenspace vectors are re-expressed in the standardized coordinates** from their stored per-module deltas (`scripts/12_restandardize_reference_vectors.py`, idempotent; cosine with the old vectors 0.97–0.99, magnitudes larger). Cosines against them from earlier versions are not comparable. Build scripts 04/06/07 use `project_shift`.
- `score_modules` coverage frame gains `aliases_used`, and a `UserWarning` names catalog genes that exact matching misses but their current symbol would recover (see Added).
- `load_attractor_direction()` / `attractor_v1.json` (the entry that was unreleased before this one cited the paper under its working title, *Cancers of all lineages converge on the placental cytotrophoblast*; the published title is *Diverse cancers converge on one region of normal cell-state space in a cancer-blind coordinate system*): the published refinement now states the paper's current result (region robust, anchor catalog-conditional; megakaryocyte a secondary target that falls with transcriptome complexity; no baseline gating) and marks the 8-cell cluster as the superseded Paper 4.2 estimate.

### Changed — documentation and examples
- README rewritten around the paper; `docs/SCOPE.md`, `docs/QUICKSTART.md` and `docs/EIGENSPACE_EVOLUTION.md` updated (the axes as the paper describes them, the standardized projection, what is retracted).
- `CITATION.cff`: preferred citation is now the final paper (Zenodo 10.5281/zenodo.20542130, a preprint) instead of the first draft; software version and DOIs unchanged.
- `CONTRIBUTING.md`: the criteria for a new module are the paper's four operational criteria documented as in Table S9b (an *unclear* entry is acceptable, the named feedback element is checked against the scored genes, nothing is written back retrospectively); new reference vectors must be in the standardized coordinates (`project_shift`); a section on statements about the convergence result.
- `ROADMAP.md`: the Paper 3 alignment is the next release; the follow-ups the paper leaves open are listed; later releases renumbered (proposal for the maintainer).
- `examples/tutorial.ipynb`: Workflow 1 describes the axes as the paper does and no longer attributes a convergence axis to PC4; Workflow 3 compares a *shift* (`project_shift`) with the reference directions; pointers to the convergence recipe. `examples/recipe_cancer_attractor.ipynb`: marked as the legacy Paper 4.2 test, eigenspace mode via `project_shift`. Builders: `scripts/08_…`, `scripts/10_…`.
- `pyproject.toml`: description and keywords.

### Added
- **`pct.convergence`** — the paper's assay: `assign_cells` / `assign_to_reference` (QC with ≥ 1,000 detected genes, `normalize_total(1e4)`+`log1p`, 44-module scoring, standardization, top-1 nearest of 152 normal cell types by cosine, erythrocyte and platelet quarantined before assignment, shape-only variant), `core_fraction` / `region_occupancy` / `enrichment`, `within_patient_occupancy` (paired Wilcoxon + sign test with the paper's conventions, ≥ 20 cells per arm, optional complexity matching), `recurrence_ranking`, `core_signature` / `module_contributions`, `prepare_expression`.
- **Bundled paper reference** (`pct.load_convergence_reference()`): definitions of the anchor / core / four-cell neighbourhood / 14-cell region, QC rules, cohort-level within-patient results (`published_within_patient()`), the retractions and limits, and Supporting-Information Tables S1, S2, S4, S6, S7, S8, extracted by `scripts/14_extract_paper_tables.py` with the totals the paper states asserted (448,399 QC cells, 574,451 assignments, contributions summing to 0.711, occupancy ≈ 100%).
- **Pooling, occupant comparison, contraction** — `pct.pooled_within_patient` (Note S-D: the per-patient malignant-minus-normal difference with a dataset random intercept, REML, Wald z and the paper's one-sided P; agrees with `statsmodels`), `pct.occupant_contrast` (Fig. 6 / S18: malignant versus matched-normal occupants of the core module by module — paired Wilcoxon over patients, detection-depth balance, Spearman against the core signature; the manuscript does not state how the sixteen modules were selected, so the table carries P and a BH q and the threshold is the caller's), `pct.contraction` and `pct.convergence.contraction_null` (Fig. 2 / Note S-G: mean pairwise cosine distance among the cancers' centroids and among their origins, the matched-averaging null in generic form with the paper's permutation P). `assign_cells` also returns `n_genes`. Table S8's statement that three of the 300 origin pairs share an origin is reproduced.
- **Reference calibrations** (`pct.load_calibration`, `pct.Calibration`, `pct.hpa_reference_scores`, `pct.read_hpa_single_cell_type`): `pseudobulk` bundled; `single_cell_scaled` buildable from the HPA table with `scripts/13_build_single_cell_calibration.py`, which publishes it only if the megakaryocyte / cytotrophoblast z-profiles reproduce Table S7. `pct.project_shift` projects module differences.
- **Gene-symbol handling (Note S-Q):** `pct.aliases_of`, `score_modules(resolve_aliases=True)` for the renamed MB21D1/TMEM173/ARNTL (→ CGAS/STING1/BMAL1) and a few renamed activity-list genes. Off by default because the bundled reference was built by exact matching (cGAS-STING on 3 genes, Circadian on 4); the paper lists re-derivation with the corrected catalog as not completed.
- **`pct.module_criteria()`** — Table S9b: sensor, TF output, feedback element, receptor field and confidence for all 44 modules (as printed; abbreviated), with the paper's counts asserted (38 feedback elements named, 6 unclear, receptor field blank for 16, ten named feedback elements outside the scored genes).
- `scripts/15_paper_reproduction_report.py` — which printed numbers a calibration reproduces; rerun it with `--calibration single_cell_scaled` after building that calibration.
- Documentation: [`docs/CONVERGENCE.md`](docs/CONVERGENCE.md), [`docs/PAPER_MAP.md`](docs/PAPER_MAP.md); a new section in [`docs/SCOPE.md`](docs/SCOPE.md) with the paper's validated, retracted, suggestive and unsupported claims; `examples/recipe_convergence_occupancy.ipynb`.
- README figure rebuilt as the paper's objective plane (cosines to the two core cells) with the published convergence ranking, light and dark, plus its data table (`scripts/09_build_eigenspace_figure.py`).
- Tests: 131 new (205 total) — statements of the paper the bundled reference reproduces (Kaiser 9; 73.8% of variance; Fig. S1 stability tiers; NPAS4 leave-one-out; Table S9 gene lists; twelve shared genes over 17 modules), the standardized-projection invariants, the paper's exact statistical conventions, the full per-cell pipeline on cells synthesized from the bundled reference, every transcribed table against the totals the paper states, and documentation integrity (relative links resolve, notebooks are valid and carry no error output, documented scripts exist).

### Fixed
- Catalog: Hippo's feedback entry (NKD1, NKD2, BIRC2, BIRC3, BIRC5) withdrawn as the paper states — NKD1/NKD2 are Wnt antagonists and the BIRCs YAP targets, not feedback inhibitors.
- `docs/EIGENSPACE_EVOLUTION.md` attributed the drop from 12 to 9 components partly to NPAS4; the paper's leave-one-out (Note S2) shows a 43-module rebuild reproduces all nine components (cosine ≥ 0.95).
- `examples/tutorial.ipynb`: the synthetic-shift demo drew its noise before seeding the generator, so its printed cosine changed on every run; it is now seeded (output regenerated: +0.931).
- Stale text: README "73 tests" (74 existed) and "This is v0.2.0"; `eigenspace` docstrings still describing 12 PCs; `reference/__init__.py` saying the disease and aging vectors were pending recomputation (done in 0.2.0); `project`'s default `n_pcs=12` against a 9-PC reference.

### Corrected statements (the paper retracts or qualifies them)
- "Malignancy deepens the state ~6×" → the effect is reported in percentage points (about +10; the fold-change is inflated by a small normal denominator).
- "kidney + liver 77/86" → 77 of 86 is the four-dataset pool (colorectal ×2, kidney, liver).
- "Absent in high-baseline lung" → **retracted**: a cell-selection artifact; clean malignant lung cells deepen (10 of 10 patients). No baseline-gating claim is made.
- "25 cancers converge on the cytotrophoblast" → the convergence *region* is the robust finding; the cytotrophoblast is the anchor under this catalog (independent catalog reconstructions name the gastric progenitor).
- "Megakaryocyte, a lineage-restricted secondary" → carried as a secondary target because it falls with transcriptome complexity like the quarantined erythrocyte, not because it is blood-specific.
- "PC1 = perception breadth; PC4 = cancer convergence axis" → the paper describes PC1 as inflammatory → terminal germline, PC2 as quiescent → proliferative, PC3 as signaling-dominated immune → metabolic epithelia, and convergence as a high-dimensional nearest-neighbour property, not an axis or a 2-D location.
- Convergence is "not convergent evolution", and a static convergent endpoint is not a proven dynamical attractor.

### Removed
- `examples/figures/eigenspace_pc1_pc4.png` (replaced by `convergence_region.png` / `convergence_region_dark.png`).

### Not done (and why)
- The `single_cell_scaled` calibration is not bundled: it derives from the HPA single-cell-type table, which is not redistributed (and was not reachable while this was prepared). Several printed numbers of the paper (Fig. 1B variance labels, Fig. S8, Table S7, the four-cell neighbourhood cosine) are computed on it and are therefore not reproduced digit for digit by the bundled pseudobulk matrix — see `docs/PAPER_MAP.md`.
- Not implemented: the specification curve of the shape contraction (the manuscript does not define its variants exactly enough to code them), the paper's seven-tissue normal pool for the matched-averaging null (the generic null needs your normal single-cell data), and the Perturb-seq / immune / reprogramming analyses (their results are recorded in the convergence reference; scripts E3–E44 are in the paper's Zenodo record).
- Notebook outputs: `examples/tutorial.ipynb` and `examples/recipe_cancer_attractor.ipynb` download PBMC3K, which could not be fetched while this was prepared. The cells whose code or output changes with the standardized projection (the tutorial's projection cell and reference comparison; the recipe's eigenspace cell) are stored without output; the others keep the outputs of the previous run, which the changes do not affect. Regenerate with `scripts/08_build_tutorial_notebook.py` / `scripts/10_build_recipe_attractor.py` followed by `jupyter nbconvert --execute --inplace`. `examples/recipe_convergence_occupancy.ipynb` (`scripts/16_build_recipe_convergence.py`) needs no download and ships fully executed.

## [0.2.2] — 2026-05-10

Discoverability + community-readiness patch. No code or data changes.

### Added

- **`CONTRIBUTING.md`** — explicit PR process with bar varying by contribution type. Bug fixes and docs casual; new canonical modules require the four perceptual criteria + at least one published a-priori prediction reproduced as a regression test in `tests/test_integration.py`. Operations falsified by Paper 4.1 (drug-disease cosine matching, TF-autoregulation as 1st-class layer, equilibrium framings, FDA-prediction by clean-signature score, snapshot-pulsatile classification) are explicitly listed as rejected — reintroducing them would invite users to do the analyses the data already refused.
- **README hero figure** (`examples/figures/eigenspace_pc1_pc4.png`, replaced in the Paper 3 alignment release) — 154 HPA cell types projected to PC1 × PC4 with the 8-cell cancer attractor cluster marked. Demonstrates the central Paper 4.2 finding visually before any prose.
- **Two focused use-case recipes** (5-10 cells each, single workflow):
  - [`examples/recipe_cancer_attractor.ipynb`](examples/recipe_cancer_attractor.ipynb) — does my (tumor − normal) shift align with the cancer-transformation direction?
  - [`examples/recipe_experimental_design.ipynb`](examples/recipe_experimental_design.ipynb) — given a cell type and a pathway, will it ramp up under stimulus, or is it saturated?
- Build scripts: `scripts/09_build_eigenspace_figure.py`, `scripts/10_build_recipe_attractor.py`, `scripts/11_build_recipe_experimental_design.py`.

### Changed

- README adds "What you get — the eigenspace at a glance" section embedding the hero figure with caption; "Try it in 5 minutes" enumerates tutorial + recipes; "Contributing" section added linking to `CONTRIBUTING.md`.

## [0.2.1] — 2026-05-10

Documentation accessibility patch. No code or data changes.

### Changed

- README rewritten for newcomer-perspective accessibility. Removed implicit assumptions that the reader knows what "perceptome", "Paper 3", "Paper 4.1-4.8", "Sun 2021", "Block 5 v1.2", or other internal-research notation mean. Restructured around three plain questions ("where does this cell sit", "what can it do", "did the measurement reach biology"). Internal references like "Paper 4.5 v1.2 amendment", "Block 5 v1.2", and "P3 PASS" replaced with descriptive language plus pointers to `docs/SCOPE.md` for the audit trail.
- Added "What problem does this solve?" section with concrete questions other tools don't directly answer.
- Added "When to use it / when not" with explicit good-fit and not-good-fit examples.
- Elevated tutorial notebook to "Try it in 5 minutes" position above technical content.
- Tests badge updated 42 → 73.

### Added

- [`examples/tutorial.ipynb`](examples/tutorial.ipynb) — executable end-to-end tutorial on PBMC3K. Auto-downloads, runs in <2 minutes, ships with output cells populated. Demonstrates all 4 workflows (geometry, perceptivity, reference comparison, drug screen) plus validity scorecard. (Originally landed after v0.2.0 tag; bundled formally in v0.2.1.)
- `scripts/08_build_tutorial_notebook.py` — re-execute on every release with `jupyter nbconvert --execute --inplace`.
- Workflow 2 caveat documented in tutorial: `A_threshold=2.0` default is calibrated for HPA pseudobulk scale; single-cell `mean_raw` scoring needs lower threshold (`A_threshold=0.3` recommended for tutorial-scale data). `predict_engagement` on the bundled HPA reference uses correct default thresholds.

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
