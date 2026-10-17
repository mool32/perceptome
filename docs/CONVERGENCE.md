# Convergence workflow — placing a cohort in the normal-cell space

This page is the user guide for `pct.convergence`, the implementation of the assay of Spiro (2026), *Diverse cancers converge on one region of normal cell-state space in a cancer-blind coordinate system* ([Zenodo 10.5281/zenodo.20542130](https://doi.org/10.5281/zenodo.20542130)). Run it on a cohort with paired tumour and adjacent-normal single-cell data; the paper released it as the general, falsifiable assay of its claim (**malignant epithelium occupies a region of normal cell-state space more often than the same patient's normal epithelium**).

For what the paper does and does not support, read [`SCOPE.md`](SCOPE.md) first. For the map from the paper to the code, see [`PAPER_MAP.md`](PAPER_MAP.md).

## What the assay does

1. **Quality control.** Gene symbols uppercased and de-duplicated; cells with fewer than 1,000 detected genes dropped (low-complexity cells carry degenerate module scores that map to low-complexity targets); `normalize_total(1e4)` + `log1p` if the matrix still holds counts.
2. **Score** the 44 regulatory modules (`mean_raw`: mean log-normalized expression of the module's genes).
3. **Standardize and project.** Each module is z-scored against a reference calibration and the result is projected through the fixed 44 × 9 loading matrix. Cancer contributes nothing to the dimensions or the geometry.
4. **Assign** each cell to the nearest of the **152 candidate normal cell types** by cosine. The erythrocyte and the platelet — post-cellular, low-complexity bodies — are removed from the pool *before* this step; the cells that would have gone to them are reassigned to their best surviving neighbour.
5. **Occupancy.** The fraction of a cell population whose nearest normal cell type is in the **core** (cytotrophoblast + megakaryocyte), or in a wider region; enrichment over the `k/152` chance baseline.
6. **Within-patient test.** Per patient, the core fraction of the malignant arm against the matched-normal arm (≥ 20 cells in each), compared across patients by a paired Wilcoxon test and a sign test. *The unit of replication is the patient, not the cell.*

```python
import scanpy as sc
import perceptome as pct

adata = sc.read_h5ad("cohort.h5ad")        # malignant + matched-normal epithelium, HGNC symbols in var_names
res = pct.convergence.assign_cells(adata, calibration="single_cell_scaled")
a = res["assignment"]                      # per cell: nearest, cosine, first_choice, displaced

pct.region_occupancy(a)                    # anchor / core / neighbourhood / region: % of cells, fold over chance

test = pct.within_patient_occupancy(
    a, adata.obs["patient"], adata.obs["arm"], malignant="tumour", normal="normal")
test["per_patient"]                        # malignant / normal fraction, difference, cell counts per patient
test["n_deeper"], test["n_patients"]       # e.g. 43, 46 in the paper's colorectal cohort
test["wilcoxon_two_sided"], test["sign_test_one_sided"]
```

Compare with the paper's cohort-level results:

```python
pct.convergence.published_within_patient()     # colorectal 18.9% vs 3.2%, 43/46; pooled 77/86; lung 10/10; prostate 7/15; ...
```

### Conventions that match the paper

- Paired Wilcoxon is **two-sided** and the sign test **one-sided** in the paper (10 of 10 patients: P = 2/1024 = 0.002; 9 of 10: sign test P = 0.011; 15 of 22: 0.07). Both sidednesses are returned.
- The default region is the **two-cell core**. The referee-control notes use a three-cell core (`region="three_cell_core"`, chance 3/152 ≈ 2%); the boundary-independence check uses `anchor`, `core`, `neighbourhood` (4 cells) and `region` (14 cells).
- `shape_only=True` (in `assign_cells` / `assign_to_reference`) subtracts each cell's own mean over the 44 modules, in reference and cells alike, **before** the reference standardization, and recomputes the reference mean and SD on the centred reference. The paper's colorectal +10.5 points, lung +12.1 and gastric +3.4 are shape-only values. (Note S-A's 7–11% comes from centring applied *after* z-standardization, which removes more than the magnitude component; the paper keeps the two variants apart and does not interchange them.)
- `within_patient_occupancy(..., n_genes=..., match_complexity=True)` resamples each patient's malignant cells to the normal arm's detected-gene distribution — the control that retains about 80% of the colorectal effect.

## Choosing a calibration

The module scores of a cell and of a reference cell type are only comparable if they sit on the same scale. The paper uses two calibrations of the same HPA reference:

| Calibration | Reference scored as | Use for | Shipped |
|---|---|---|---|
| `pseudobulk` | mean log1p(nCPM) on population-averaged HPA expression | cell-type-level data on the HPA scale; what the eigenspace is built on | yes |
| `single_cell_scaled` | each HPA cell type scaled to 1e4 total, then log1p — the preprocessing of the per-cell data | **per-cell data** | no (build it) |

The single-cell definition above is this package's reading of Note S1 ("the same normalize_total(1e4)+log1p preprocessing applied to the tumor data"): `pct.hpa_reference_scores(..., scaling="single_cell_scaled")`. The Table S7 check in the build script is what confirms it. Under the pseudobulk calibration **tumour cells collapse onto low-complexity targets** (core fraction 7% instead of 17%, centroid distances degenerate), and so do independent normal single cells — 99.6% of 16,312 skin, kidney and breast cells land on erythroid targets (Notes S-B, S-H). That is a property of the calibration, not of cancer, which is why it is not a calibration chosen to favour the result. `pct.project` warns when per-cell-scale scores meet the pseudobulk calibration, and `assign_to_reference` does the same.

**Building the single-cell calibration.** The HPA single-cell-type table is not redistributed with this package (and could not be fetched while this release was prepared):

```bash
# download rna_single_cell_type.tsv from https://www.proteinatlas.org/download/ (the release used for the paper)
python scripts/13_build_single_cell_calibration.py rna_single_cell_type.tsv
```

The script rebuilds the pseudobulk reference from the table and compares it with the bundled matrix (an independent check that the table is read correctly and is the same release), builds the single-cell reference, and checks the result against the paper: the z-scores of the megakaryocyte and the cytotrophoblast over the 44 modules must reproduce Table S7. Only then is it saved under the name `pct.load_calibration("single_cell_scaled")` loads; otherwise it is written as `*.UNVERIFIED.json` and the script exits with status 2. Afterwards `python scripts/15_paper_reproduction_report.py --calibration single_cell_scaled` shows which printed numbers of the paper it reproduces.

## Reading the result

What the paper establishes, in its own terms:

- **The region is the robust finding; the anchor is catalog-conditional.** Three independent reconstructions of the catalog from the same four criteria reproduce the convergence region (five normal cell types are in the top ten of all four catalogs) but name the gastric progenitor, not the cytotrophoblast, inside it. Quote the cytotrophoblast as "the anchor under this catalog".
- **The within-patient increase is the most defensible claim**, and the most portable: KEGG, Reactome, Hallmark and PROGENy collections put through the same pipeline reproduce it. What this catalog adds is the cross-lineage generality of the target.
- **Convergence is not convergent evolution.** The data are snapshots; the test shows a convergent endpoint across origins, not a path, and not a proven attractor.
- **Most of the raw contraction is magnitude.** With the overall-magnitude component removed a smaller, real convergence in module shape remains — report `shape_only` alongside raw values.

Pitfalls the paper documents:

- **Malignant-cell purity.** A tumour arm that is mostly normal-adjacent epithelium is diluted to the normal baseline (the Census lung "tumour" arm was ~9% annotated-malignant, which produced the retracted "no deepening in lung"). Use author-annotated malignant compartments; marker-based selection biases the other way.
- **The right comparator.** Colorectal cancer arises from crypt-base progenitors, not from the tissue at large; against differentiated surface epithelium the effect is overstated (+19.1 points vs +14.8 against the crypt base).
- **Do not rank cancers by convergence fraction.** The ordering is sensitive to malignant-label quality (an author-validated relabeling moved one colon dataset from 13% to 36%) and to residual platform effects.
- **Report points, not folds.** The fold-change is inflated by a small normal denominator; on the additive scale the increase is roughly constant (~+10 points) across tissues.
- **Pooling.** Pooling patients across datasets can manufacture a Simpson's-paradox artifact; the paper uses a mixed model with a dataset random intercept (+13.0 points over four datasets). `pct.pooled_within_patient` fits it (see below).
- **Gene symbols.** Three catalog genes have been renamed (MB21D1, TMEM173, ARNTL → CGAS, STING1, BMAL1). With current symbols cGAS-STING is scored on three genes and Circadian on four, matching the bundled reference; `resolve_aliases=True` recovers all five but then differs from the reference until it is re-derived.

## Pooling cohorts

Several datasets of the same question are pooled with the paper's mixed model (Note S-D): the per-patient malignant-minus-normal difference, an intercept, and a random intercept per dataset (REML, Wald z, one-sided P for a deepening):

```python
tests = {name: pct.within_patient_occupancy(a, patient, arm, malignant="tumour", normal="normal")
         for name, (a, patient, arm) in cohorts.items()}
pooled = pct.pooled_within_patient({name: t["per_patient"] for name, t in tests.items()})
pooled["estimate_points"], pooled["se_points"], pooled["z"], pooled["p_one_sided"]    # paper: +13.0, 1.87, 6.95, 1.9e-12
pooled["per_dataset"]                                                                  # patients, patients deeper, mean difference
```

With few datasets the between-dataset variance is estimated coarsely (the paper's four are two colorectal cohorts, kidney and liver), and the one-sided P does not propagate the uncertainty of that estimate; say so when you report it.

## What malignancy adds inside the state (Fig. 6)

Among the cells that already occupy the core, how does a malignant cell differ from a normal cell of the same patient? `occupant_contrast` compares the two occupant populations module by module — per patient (≥ 10 occupants in each arm), a paired Wilcoxon test across patients, and a detection-depth balance (the malignant occupants resampled to the normal occupants' detected-gene distribution, as in Note S-F):

```python
res = pct.convergence.assign_cells(adata, calibration="single_cell_scaled")
out = pct.occupant_contrast(res["scores"], res["assignment"], adata.obs["patient"], adata.obs["arm"],
                            malignant="tumour", normal="normal", n_genes=res["n_genes"])
out["table"]                       # per module: raw and depth-balanced difference, patients higher, P, BH q
out["rho_core_signature"]          # Spearman over the modules against the core's defining signature (paper: +0.36)
out["depth"]                       # detected genes of the occupants, raw and balanced (paper: 4,351 vs 2,395 -> 2,640 vs 2,496)
```

The paper reports this for colorectal only ("awaits replication") and names sixteen modules that separate the occupants; it does not state the rule that selected them, so the table gives P and a Benjamini-Hochberg q (an addition of this package) and the threshold is yours. The paper's reading — malignancy superimposes stress, growth and metabolic programs on the shared state rather than moving cells further along the axis that defines it — is a statement about that colorectal comparison, not a property this function establishes.

## Contraction

The paper's cross-origin contraction pairs each cancer's malignant centroid (the mean of its quality-controlled cells in the nine-dimensional space) with its normal cell of origin (Table S8) and compares the mean pairwise cosine distance among the origins (about 0.99: near-orthogonal states) with that among the malignant centroids (about 0.16):

```python
res = pct.convergence.assign_cells(adata, calibration="single_cell_scaled")
origins = pct.load_convergence_reference()["tables"]["S8"]["cell_of_origin"]      # cancer -> normal cell of origin
c = pct.contraction(res["coordinates"], adata.obs["cancer"], origins=origins, calibration="single_cell_scaled")
c["mean_distance"], c["origin_mean_distance"], c["ratio"]
```

The paper takes this headline apart, and so should you. The raw contraction is largely a matched-averaging and shared-magnitude effect: aggregates of unrelated normal cells formed by the same averaging reach about 0.13 (`pct.convergence.contraction_null`, which needs normal single-cell data of your own, projected with the same calibration, and your definition of a coherent micro-population). The comparison with random draws of 25 normal cell types (a Stayton-type test) and the fold-contraction it produced are withdrawn. What survives is a small convergence in module shape; the paper reports it as a specification curve because the result depends on how magnitude-removed scores are rescaled (two of three standardizations support it, one does not). Run `contraction` on `shape_only=True` coordinates and on the other variants, and report them side by side.

## What is not implemented here

The 36-way specification curve (the manuscript does not define its within-cell-centred and per-cell-rank variants or the third standardization exactly enough to code them), the Perturb-seq, immune and reprogramming analyses, and the independent-catalog reconstructions. The paper's seven-tissue normal pool for the matched-averaging null is not part of the package either. Results of the analyses are recorded in `pct.load_convergence_reference()`; the analysis scripts (E3–E44) are in the Zenodo record of the paper.
