# Paper ↔ repository map

The paper: Spiro (2026), *Diverse cancers converge on one region of normal cell-state space in a cancer-blind coordinate system* ([Zenodo 10.5281/zenodo.20542130](https://doi.org/10.5281/zenodo.20542130)), main text and Supporting Information of manuscript version v11. It is "Paper 3" of the series this repository implements; Papers 4.x (memory consolidation, drugs, capacity layer) are separate and are not covered by it.

This page says, for each part of the paper, **where it lives in the code and how far it can be checked**:

| Status | Meaning |
|---|---|
| **reproduced** | a test recomputes the printed number from the bundled data |
| **implemented** | code exists and is tested on synthetic data built from the bundled reference |
| **bundled** | the paper's own numbers are shipped as data (checked against totals the paper states) |
| **external** | needs data that are not redistributed here (HPA single-cell-type table, cancer cohorts) |
| **not implemented** | analysis not part of this package |

## The coordinate system

| Paper | Where | Status |
|---|---|---|
| 44 modules, 206 unique genes (219 slots), Table S9 | `catalog/data/modules_v03.json` (`core_genes`); `tests/test_paper_reference.py` | **reproduced** — identical, gene by gene |
| Twelve shared genes touching 17 modules; definitional pairs (XBP1, STAT1) — Note S-L | same test | **reproduced** |
| Four operational criteria per module, Table S9b (38 feedback elements named, 6 unclear, receptor field blank for 16, ten named feedback elements outside the scored genes) | `pct.module_criteria()`; `catalog/data/module_criteria_s9b.json` | **bundled** as printed (gene lists abbreviated in the print); counts **reproduced**. Hippo's withdrawn feedback entry is removed from the catalog |
| Nine components by the Kaiser criterion, 73.8% of variance; bootstrap stability 0.63–0.99 (PC1–PC6) and 0.40–0.54 (PC7–PC9), PC8 below PC9 — Note S1, Fig. S1 | `eigenspace/data/reference_v03.json` | **reproduced** (rebuilding from the bundled 154 × 44 matrix gives the shipped reference exactly) |
| NPAS4 leave-one-out: 43-module rebuild reproduces all nine components, per-PC cosine ≥ 0.95, mean 0.99 over PC1–PC6 — Note S2 | `tests/test_paper_reference.py` | **reproduced** (min 0.951, mean 0.993) |
| PC1 runs from inflammatory to terminal germline cells | `pct.project` | **reproduced** in sign and ordering (neutrophils, monocytes, macrophages negative; oocytes, megakaryocytes positive) |
| Unrelated normal states are near-orthogonal (cosine distance ~0.99) | `pct.project` | **reproduced** qualitatively (mean pairwise cosine 0.002 over the 154 types; 0.976 distance among the 25 origins) |
| z-score against the reference, then project — Methods, Note S1 | `pct.project` (`standardize=True`), `pct.load_calibration` | **implemented**. Earlier versions skipped the z-scoring |
| Pseudobulk and single-cell-scaled calibrations — Notes S1, S-B, S-H | `pct.eigenspace.calibration`; `scripts/13_build_single_cell_calibration.py` | pseudobulk **bundled**; single-cell-scaled **external** (build it; the script publishes it only if Table S7 is reproduced) |

## The convergence assay

| Paper | Where | Status |
|---|---|---|
| Cancer projection: ≥1,000 genes/cell, `normalize_total(1e4)`+`log1p`, 44 modules, standardize, project, top-1 nearest centroid by cosine | `pct.convergence.assign_cells` / `assign_to_reference` / `prepare_expression` | **implemented** (cells synthesized from the reference return to their own cell type) |
| Erythrocyte and platelet quarantined before assignment; 152 candidates — Notes S3, S-K | same | **implemented** |
| Shape-only transform, defined and ordered — Note S-K (iv) | `shape_only=True` | **implemented** |
| Two-cell core, three-cell core, four-cell neighbourhood, 14-cell region; enrichment over `k/152` — Fig. 3F | `pct.convergence.region_members`, `region_occupancy`, `enrichment` | core and neighbourhood **reproduced** from Table S2 (23.9%, 41.9%; 18.2×, 15.9×; anchor 11.1%, 16.9×); the 14-cell region within 1 point of the printed 76.9% (see *Observations* below) |
| Cross-cancer recurrence, depth, occupancy — Table S2, Fig. 3B–C | `pct.recurrence_ranking` (your cohorts); Table S2 bundled | **implemented** / **bundled** |
| Per-module contributions of a cancer to the core signature — Table S4, Fig. 4 | `pct.convergence.core_signature`, `module_contributions` | **implemented**; Tables S4 and S7 **bundled** (S4 sums to the stated 0.711) |
| Within-patient occupancy test (paired Wilcoxon + sign test, ≥20 cells per arm) — Fig. 5 | `pct.within_patient_occupancy` | **implemented**; the paper's exact statistical conventions are **reproduced** (two-sided Wilcoxon 2/1024 for 10 of 10; one-sided sign tests 0.011 for 9/10 and 0.07 for 15/22) |
| Complexity-matched resampling — Note S-F | `match_complexity=True` | **implemented** |
| Cohort-level results: Pelka colorectal, HTA colorectal, kidney, liver, pooled, gastric, lung, prostate ×2, IPF, breast, adenoma, DCIS, Barrett | `pct.convergence.published_within_patient()` | **bundled** (internal consistency tested: pooled 77/86 = 43+23+6+5 over 46+25+7+8) |
| The 25 cohorts, cell counts, labels — Tables S1, S8; complexity of assigned cells — Table S6 | `pct.load_convergence_reference()["tables"]` | **bundled** (448,399 QC cells and 574,451 assignments reproduced) |
| Retractions and limits | `pct.load_convergence_reference()["retracted"]`, `["limits"]`; [`SCOPE.md`](SCOPE.md) | **bundled** |
| Mixed-effects model over datasets, dataset random intercept — Note S-D | `pct.pooled_within_patient` | **implemented** (intercept plus random intercept, REML, Wald z with the paper's one-sided P; agrees with `statsmodels`); the paper's +13.0 points over 86 patients needs per-patient data that are not published, so it is **not reproduced** |
| Malignant versus matched-normal occupants of the state, module by module — Fig. 6, Fig. S18 | `pct.occupant_contrast` | **implemented** from the manuscript's description (the patient as unit, ≥ 10 core-assigned cells per arm, paired Wilcoxon, detection-depth balance by nearest-gene-count matching, Spearman against the core signature). The text does not give the rule that selects the sixteen modules, so the table carries P and a Benjamini-Hochberg q (this package's addition) and the threshold is the caller's; colorectal only in the paper ("awaits replication"); the finding is recorded in [`SCOPE.md`](SCOPE.md); the paper's modules and ρ = +0.36 are **not reproduced** (no cohort data) |
| Cross-origin contraction: mean pairwise cosine distance among the cancers' malignant centroids and among their origins; "three of the 300 origin pairs share an origin" — Fig. 2, Methods, Table S8 | `pct.contraction` | **implemented**; the origin pairing is **reproduced** (300 pairs, three shared: B-ALL/DLBCL, germ-cell/seminoma, osteosarcoma/Ewing) and the origins' mean distance is 0.976 on the bundled calibration (paper 0.99); the malignant 0.16 needs the cohorts |
| Matched-averaging null of the contraction — Note S-G | `pct.convergence.contraction_null` | **implemented in generic form** (aggregates of the same sizes from coherent normal micro-populations, the paper's permutation P); it needs normal single-cell data and the caller's definition of a coherent micro-population — the paper's seven-tissue pool and its definition are not in the manuscript — so its 0.129 is **not reproduced** |
| Specification curve of the shape contraction — Note S-G, Fig. S19 (36 crossings: calibration × magnitude handling × core definition × erythroid handling) | — | **not implemented**: the manuscript does not define the within-cell-centred and per-cell-rank variants or the third standardization exactly enough to code them. `contraction` and `contraction_null` take any coordinates, so each variant can be run and reported side by side |
| Catalog-dependence reconstructions (Note S-O), Perturb-seq (S-J), immune (S-J), reprogramming (Fig. S28), time course (Fig. S26) | results in `convergence_v1.json` | **not implemented**; scripts E3–E44 are in the Zenodo record |

## Gene symbols (Note S-Q)

Auditing the scored panel found catalog genes that matched neither HPA nor most cohorts: three renamed symbols (MB21D1→CGAS and TMEM173→STING1 in cGAS-STING, ARNTL→BMAL1 in Circadian) and NKX3-1, lost by a loader that stripped `-<digits>` suffixes. Here: `pct.aliases_of`, `score_modules(resolve_aliases=...)`, a warning that names recoverable genes, and a `prepare_expression` that never strips suffixes. The paper states that the corrected catalog **has not been re-derived end to end**; neither has the bundled reference (cGAS-STING is scored on three genes and Circadian on four), which is why `resolve_aliases` defaults to off.

## What the bundled data do not reproduce

The following printed numbers are computed on the **single-cell-scaled** calibration (Fig. 1B variance labels, Fig. S8, Table S7, the neighbourhood cosine) or on a rebuilt matrix; the shipped pseudobulk reference gives different values. `python scripts/15_paper_reproduction_report.py` prints this table for any calibration; rerun it with `--calibration single_cell_scaled` after building that calibration — the rows should then turn to "match".

| Statement (paper) | Paper | Bundled pseudobulk reference |
|---|---|---|
| PC1 / PC2 share of variance (Fig. 1B) | 31% / 21% | 29.2% / 13.7% (the nine-component total, 73.8%, does match) |
| Module pairs with \|r\| > 0.7 of 946 (Fig. S8) | 6 | 20 |
| Mean \|r\| over the 946 pairs | 0.21 | 0.27 |
| The six pairs | ATF6–IRE1, JAK-STAT–Type I IFN, ERK–NF-κB, PERK–p53, NF-κB–PERK, BMP–Hippo | 5 of the 6 (JAK-STAT–Type I IFN is below 0.7) |
| Table S7 z-scores of the two core cells | printed (2 decimals) | Pearson r 0.71 (megakaryocyte) / 0.78 (cytotrophoblast); not equal |
| Mean cosine distance among the 25 origins | 0.99 | 0.976 |
| Mutual cosine within the four-cell neighbourhood | 0.77 | 0.52 |
| Manifold medoid | alveolar type 1 | respiratory ciliated cells (Euclidean, 9-D) |

## Observations for the paper's author

Found while aligning the repository; none changes the code.

1. **Fig. 1B variance labels vs the Methods.** The nine-component total (73.8%), the Kaiser count and the Fig. S1 stability tiers match the shipped `reference_v03` exactly; the PC1 / PC2 labels (31% / 21%) do not (29.2% / 13.7%). Presumably they come from a PCA of the single-cell-scaled matrix.
2. **Note S-N (iii): the rebuilt k = 44 space differs from the shipped reference (cell-to-cell r = 0.964).** Rebuilding from the bundled 154 × 44 matrix reproduces `reference_v03` to machine precision (eigenvalues to 2e-15), so the matrix behind the paper's rebuild is a different one.
3. **Fig. 3F vs Table S2.** The 14-cell region prints 76.9% (8.4×). The twelve cells with recurrence ≥ 4 sum to 73.95% in Table S2; places 13–14 are a five-way tie at recurrence 2 (enteric stem cells, erythrocyte progenitors, Paneth cells, syncytiotrophoblasts, undifferentiated spermatogonia), and depending on which two are taken the 14-cell sum is 76.7–77.9%. The printed 76.9% is attainable, but the table does not say which pair was used. `region_members("region")` breaks the tie by median then mean depth (77.6%).
4. **"Table 1"** is cited (Results, Fig. 3F; Note S-K) but is not in the manuscript.
5. **Note S-Q lists three genes lost "for two reasons"** (NKX3-1; MB21D1/TMEM173; ARNTL) — four symbols in three modules.
6. **Table S9b note: the named feedback element is absent from the scored genes "in 11 of the 38 modules that name one".** The eleven listed include Hippo, whose feedback entry the same note withdraws (it is among the six "unclear", not the 38); ten of the 38 remain, which is what `pct.module_criteria()` gives.

## The pre-paper repository

Statements that earlier versions of this repository made and the paper does not support, corrected in this release: PC1 as "perception breadth" and PC4 as a "cancer convergence axis" (the paper describes the axes as inflammatory → germline, quiescent → proliferative, immune → metabolic epithelia, and convergence as a high-dimensional nearest-neighbour property); "25 cancers converge on the cytotrophoblast" without the catalog-conditional qualification; the "~6×" fold-change; "kidney + liver 77/86" (77 of 86 is the four-dataset pool); "absent in high-baseline lung" (retracted); the megakaryocyte as a "lineage-restricted secondary" (the paper rejects it because it falls with transcriptome complexity, not because it is blood-specific).
