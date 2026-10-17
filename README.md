# perceptome

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python: ≥3.8](https://img.shields.io/badge/Python-≥3.8-blue.svg)](https://www.python.org)
[![Tests: 205](https://img.shields.io/badge/tests-205%20passing-green.svg)](#tests)
[![Software DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.20113468.svg)](https://doi.org/10.5281/zenodo.20113468)
[![Paper DOI](https://img.shields.io/badge/paper-10.5281%2Fzenodo.20542130-blue.svg)](https://doi.org/10.5281/zenodo.20542130)

A Python toolkit for the **cancer-blind coordinate system of cell state** described in Spiro (2026), [*Diverse cancers converge on one region of normal cell-state space in a cancer-blind coordinate system*](https://doi.org/10.5281/zenodo.20542130): **44 regulatory modules** — major signaling and stress-response circuits (NF-κB, mTOR, UPR, p53, nuclear receptors, …) — scored in single-cell RNA-seq data and placed in a **fixed nine-dimensional reference space built from 154 normal human cell types** (Human Protein Atlas). Cancer contributes neither to the choice of dimensions nor to the geometry, so where cancers fall in that space is an out-of-sample observation, and the analysis is falsifiable: preserved tissue identity would leave cancers dispersed near their origins.

*Terminology.* The paper calls the 44 units *regulatory modules*; earlier versions of this repository (and the substrate / drug papers it also implements) call the same objects *perceptual modules* and the framework a *perceptual system*. The API keeps its names (`pct.score_modules`, catalog categories such as `A_exteroceptive`); the prose follows the paper.

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="examples/figures/convergence_region_dark.png">
  <img alt="Two panels. A: the 154 normal human cell types of the bundled reference on the plane of proximity to the cytotrophoblast and the megakaryocyte; cell types are coloured by the share of malignant cells from 25 cancers that land nearest to them, with the cytotrophoblast (anchor) and the megakaryocyte (secondary) starred at the high-proximity edge. B: the twelve normal cell types receiving most of the convergence, led by the megakaryocyte (12.8%) and the cytotrophoblast (11.1%) against 0.66% expected by chance." src="examples/figures/convergence_region.png">
</picture>

*Panel A places the 154 normal cell types on the paper's objective plane (Fig. S6A): cosines to the two core cells in the standardized nine-dimensional space of the bundled reference. Colour and panel B are the paper's published results (Table S2), not recomputed here — the package ships no cancer data. [`convergence_region_data.csv`](examples/figures/convergence_region_data.csv) holds the plotted values; [`scripts/09_build_eigenspace_figure.py`](scripts/09_build_eigenspace_figure.py) regenerates the figure.*

For each cell or cell type, perceptome answers three questions other tools don't directly answer:

1. **Where does this cell sit among normal cell states — and do malignant cells sit closer to one region than the same patient's normal cells?** Projection into the fixed 9-D space, assignment to the nearest of 152 normal cell types, and the paper's within-patient occupancy test (`pct.convergence`).
2. **What can this cell do under perturbation?** A capacity-floor predictor: given a cell type and a pathway, will it ramp up under stimulus or is it already saturated? Validated across several experimental paradigms (intestinal homeostasis, organoid regeneration, muscle hypertrophy).
3. **Did the measurement reach biology, or hit a confound?** A built-in scorecard with three null controls (random-gene panel, housekeeping, cell cycle) that catches common artifacts in perturbation analyses before they corrupt downstream interpretation.

The tool implements **only what the underlying studies validated** and records what they falsified or retracted — see [`docs/SCOPE.md`](docs/SCOPE.md). How every statement of the paper maps onto the code, and which printed numbers the bundled data do and do not reproduce, is in [`docs/PAPER_MAP.md`](docs/PAPER_MAP.md).

---

## What the paper found

Malignant single cells of **25 cancers** — epithelial, hematopoietic, glial, neural-crest, two sarcomas, two germ-cell tumours, neuroendocrine — were projected into the fixed space. The headline numbers, each with the caveat the paper attaches to it:

| Finding | Number | Caveat stated in the paper |
|---|---|---|
| Cancers converge on a restricted region | mean pairwise cosine distance **0.16** among the malignant states vs **0.99** among their near-orthogonal normal origins | Most of the contraction is shared transcriptional output and the averaging itself (matched-averaging null 0.13). A smaller convergence in module *shape* remains: below all 200 matched draws in raw module space, P = 0.03 under the reference scaling, absent under a third scaling |
| The nearest normal cell is the placental cytotrophoblast | a top target in **18 of 25** cancers; **11.1%** of malignant cells vs 0.66% by chance = **16.9-fold** | The two-cell core adds the megakaryocyte (23.9%, 18.2-fold), which falls with transcriptome complexity exactly as the quarantined erythrocyte does and is carried as a secondary target only. Independent reconstructions of the catalog reproduce the **region** but name the gastric progenitor instead: *the region is the robust finding; the anchor is conditional on this catalog* |
| Malignant epithelium occupies the state more than the same patient's normal epithelium | colorectal **18.9% vs 3.2%**, 43 of 46 patients, paired P < 10⁻¹²; **+10.5 points** with magnitude removed (44/46); pooled over four datasets 77 of 86 patients (+13.0 points, mixed model); held-out gastric 18/22; lung 10/10 | Present in every tissue tested **except prostate**, an open negative (its two cohorts disagree about the baseline itself). Not in fibrotic lung (IPF +3.8 points, n.s.). The paper calls this within-patient increase — not the cross-lineage convergence — its most defensible claim |
| The state is normal physiology | healthy lung epithelium 16.1%, neural tissue 7.6%, kidney 0.5% | So the unit of the result is a *physiological cell state whose occupancy malignancy increases*, not a cancer-specific destination |

**What the paper does not claim.** Convergence is a static geometric test — a convergent endpoint across origins, *not* a measured trajectory and not a proven dynamical attractor. An earlier reading of the increase as gated by tissue baseline (absent in high-baseline lung) is **retracted** (a cell-selection artifact), as is the fold-change framing ("~6×"; the paper reports points). The state is not moved by single-gene perturbation (Perturb-seq in K562), malignant cells in it show no detectable link to immune infiltration (62 colorectal tumours), so the trophoblast's immune privilege does not appear to be inherited, and the pre-invasive prediction fails in esophagus. The full list — retractions, limits and open items — is in [`docs/SCOPE.md`](docs/SCOPE.md) and in `pct.load_convergence_reference()`.

---

## Try it in 5 minutes

- 📓 [`examples/recipe_convergence_occupancy.ipynb`](examples/recipe_convergence_occupancy.ipynb) — the paper's within-patient assay end to end on a synthetic paired cohort (offline, outputs populated). **Start here for the cancer analysis.**
- 📓 [`examples/tutorial.ipynb`](examples/tutorial.ipynb) — walkthrough of the geometry, capacity, reference-comparison and drug workflows on PBMC3K (downloads the demo dataset).
- 📖 [`docs/QUICKSTART.md`](docs/QUICKSTART.md) and [`docs/CONVERGENCE.md`](docs/CONVERGENCE.md) — the same workflows in markdown, with notes on common pitfalls.
- 🍳 Focused recipes: [experimental design](examples/recipe_experimental_design.ipynb) (will a pathway ramp up in my cell type, or is it saturated?) and the [legacy cancer capacity-direction recipe](examples/recipe_cancer_attractor.ipynb).

```bash
pip install git+https://github.com/mool32/perceptome.git
```

```python
import scanpy as sc
import perceptome as pct

adata = sc.read_h5ad("cohort.h5ad")       # malignant + matched-normal epithelium; HGNC gene symbols

# paper pipeline: >=1,000 genes/cell -> 44 module scores -> standardize -> project -> nearest of 152 normal cell types
res = pct.convergence.assign_cells(adata, calibration="single_cell_scaled")

# the paper's falsifiable general assay: is the tumour arm deeper in the core than the same patient's normal arm?
test = pct.within_patient_occupancy(
    res["assignment"], adata.obs["patient"], adata.obs["arm"], malignant="tumour", normal="normal")
print(f"{test['n_deeper']} of {test['n_patients']} patients deeper; "
      f"paired Wilcoxon P = {test['wilcoxon_two_sided']:.2g}")
```

**Per-cell data need the `single_cell_scaled` calibration.** Under the bundled population-averaged (`pseudobulk`) calibration, cells collapse onto low-complexity reference types — the paper shows it for tumour cells and for independent normal single cells (99.6% onto erythroid targets). The single-cell calibration is derived from the HPA single-cell-type table, which is not redistributed here: build it once with [`scripts/13_build_single_cell_calibration.py`](scripts/13_build_single_cell_calibration.py), which publishes it only if it reproduces the paper's Table S7. See [`docs/CONVERGENCE.md`](docs/CONVERGENCE.md).

The geometry workflow on its own, and the compare-to-references step:

```python
scores = pct.score_modules(adata)["scores"]         # 44 modules per cell (mean_raw)
coords = pct.project(scores)["coordinates"]         # standardized 9-PC reference coordinates
refs = pct.compare_to_references(coords)            # cosines vs legacy cancer / disease / aging directions
```

There are three more workflows — see the tutorial.

---

## What problem does this solve?

Most single-cell analyses operate at the gene level (DEGs, marker genes, GSEA on gene sets) or at the cluster level (cell type labeling, trajectory inference). Both are useful. Neither directly answers questions like:

- *"Do my tumour cells sit in a particular region of normal cell-state space, and more so than the same patient's normal cells?"*
- *"My cells of type X — can they ramp pathway Y under stimulus, or is that pathway already at ceiling?"*
- *"My perturbation effect on module Z — is it real biology or a baseline artifact from the cells being more proliferative?"*

perceptome works one level up: **44 modules, read through their induced transcriptional output**, projected into a fixed coordinate system where any cell can be compared to any reference — plus a capacity layer that predicts the ramp room of each pathway in each cell type before you run the experiment.

## When to use it / when not

**Good fit:**
- You have scRNA-seq (or cell-type-level RNA-seq) and want a low-dimensional, interpretable view of regulatory-module activity in a frame that was built without your disease of interest
- You have paired tumour / adjacent-normal single-cell data and want to run the within-patient occupancy test, or you want to see which normal cell type each malignant cell resembles most
- You're designing or interpreting a perturbation experiment (drug, knockdown, stimulus) and want a-priori predictions of which pathways have ramp room in which cell types
- You need a built-in validity check before trusting perturbation effects

**Not a good fit:**
- Single-gene differential expression — use `scanpy.tl.rank_genes_groups` or DESeq2 directly
- Cell type identification or label transfer — use scanpy / CellTypist / Symphony. (The assignment to normal cell types here is a geometric probe of one region of the space, not an annotation method; the paper contrasts the fixed linear frame with label transfer.)
- Trajectory or dynamics — the convergence result is a static geometric one; use Slingshot / scFates / palantir for trajectories
- Drug-disease repositioning via geometric similarity — this was tested in three independent ways and falsified; the tool deliberately doesn't support it (see [`docs/SCOPE.md`](docs/SCOPE.md))
- Phosphorylation cascade readout — the framework operates on transcription, not protein phosphorylation
- Clinical decision support — research tool only, no clinical validation

---

## What's in the box

| Subpackage | Function | What it answers |
|---|---|---|
| **catalog** | `pct.list_modules()`, `pct.get_genes(module, gene_set)`, `pct.module_criteria()`, `pct.aliases_of` | What modules and gene panels are defined? What are their four operational criteria (Table S9b)? |
| **score** | `pct.score_modules`, `pct.score_readiness`, `pct.score_activity` | Module activity per cell from expression |
| **eigenspace** | `pct.project`, `pct.project_shift`, `pct.load_calibration`, `pct.rebuild` | Standardize against a reference calibration and project into the fixed 9-PC space |
| **convergence** | `pct.convergence.assign_cells`, `assign_to_reference`, `core_fraction`, `region_occupancy`, `within_patient_occupancy`, `pooled_within_patient`, `occupant_contrast`, `contraction`, `recurrence_ranking`, `core_signature`, `load_convergence_reference` | Which normal cell type is each cell nearest to; how much of a cohort sits in the core; is the tumour arm deeper than the matched normal arm (and pooled over datasets); how do malignant and normal occupants of the state differ; how tightly do cancers contract relative to their origins |
| **perceptivity** | `pct.compute_perceptivity`, `pct.predict_engagement`, `pct.capacity_floor` | Capacity, headroom, saturation per cell type × module |
| **compare** | `pct.compare_conditions`, `pct.compare_to_references`, `pct.infrastructure_regime`, `pct.divergence_score` | Disease vs control deltas, regime classification |
| **reference** | `pct.load_attractor_direction`, `pct.attractor_alignment` | Legacy cancer capacity-direction (Paper 4.2) |
| **validity** | `pct.validate_perturbation`, `pct.scorecard` | Random-200 / housekeeping / cell-cycle nulls |
| **drugs** | `pct.drug_anchors`, `pct.activity_layer_screen` | 9 validated mechanism-pathway anchors (narrow scope) |
| **network** | `pct.compute_network`, `pct.module_heterogeneity` | Module co-variation + rich-club analysis |

## Bundled reference data

- **44-module catalog** — `core_genes` are the 206 unique genes of the paper's Table S9 (219 gene slots; Calcium has 4 genes, every other module 5); activity / sensor / cascade / TF / feedback lists for the capacity layer; and the four operational criteria per module as printed in Table S9b (`pct.module_criteria()`)
- **154 × 44 cell-type-by-module reference matrix** computed from the Human Protein Atlas (mean log1p(nCPM), the *pseudobulk* calibration), with R (readiness), A (activity), C (capacity) and headroom precomputed
- **9-PC reference eigenspace** built from that matrix by correlation PCA — nine components by the Kaiser criterion (73.8% of variance), PC1–PC6 bootstrap-stable, PC7–PC9 probable — and its pseudobulk **calibration** (module means / SDs, 154 normal-cell centroids)
- **Convergence reference** (`pct.load_convergence_reference()`) — the anchor / core / neighbourhood definitions, the quality-control rules, the cohort-level within-patient results, the retractions and limits, and the paper's Supporting-Information Tables S1, S2, S4, S6, S7 and S8 as DataFrames
- **Legacy cancer capacity-direction** — the Paper 4.2 per-module Δ vector and 8-cell attractor cluster from an 11–13 cancer cohort. The paper's convergence result supersedes the 8-cell cluster as the description of where cancers land; the vector stays available for the capacity-space alignment test
- **Disease reference vectors** for rheumatoid arthritis, Alzheimer's, idiopathic pulmonary fibrosis and diabetic kidney disease, and **aging reference vectors** (inflammaging: blood old−young; collapse: bone marrow old−young) — recomputed from raw scRNA-seq and expressed in the standardized coordinates. The paper uses the same four diseases as negative controls and excludes rheumatoid arthritis because its atlas has a median of 286 detected genes per cell, below the 1,000-gene quality floor — check the complexity of the RA data behind this vector before relying on it
- **9-anchor drug reference table** with full validation chain
- **Validity panels**: housekeeping (Eisenberg & Levanon 2013), cell cycle (Tirosh 2016)

## Methodology notes

- **Standardized projection.** A cell is placed in the space by scoring the 44 modules (`mean_raw`), z-scoring each against a reference calibration, and projecting through the fixed loadings (`pct.project`). Versions ≤ 0.2.2 multiplied raw scores by the loadings, which left every cell sharing one large offset; `project(..., standardize=False)` still reproduces that, but it is not the paper's geometry.
- **Two calibrations.** `pseudobulk` (bundled; what the eigenspace is built on) is for cell-type-level, population-averaged data. `single_cell_scaled` — the reference re-derived under the same `normalize_total(1e4)` + `log1p` preprocessing as the per-cell data — is required for per-cell data. `pct.project` warns when per-cell-scale scores meet the pseudobulk calibration. Several printed numbers of the paper (Fig. 1B variance labels, Fig. S8, Table S7) are computed on the single-cell calibration and are therefore not reproduced digit for digit by the bundled matrix; `scripts/15_paper_reproduction_report.py` lists them.
- **Shape-only.** Subtracting each cell's own mean over the 44 modules (reference and cells alike) removes overall regulatory magnitude; most of the raw contraction is magnitude, and the within-patient increase survives its removal (`shape_only=True`).
- **Quality control.** Cells need ≥ 1,000 detected genes: low-complexity cells carry degenerate module scores that map to low-complexity targets. The erythrocyte and platelet are quarantined from the 154 reference cells before assignment (152 candidates). Report the unit of replication as the patient, not the cell.
- **Gene symbols.** Catalog genes are matched by exact symbol, as the bundled reference was built. Three catalog genes have been renamed (MB21D1→CGAS, TMEM173→STING1, ARNTL→BMAL1), so data annotated with current symbols score cGAS-STING on three genes and Circadian on four; `score_modules` warns and `resolve_aliases=True` recovers them. The paper lists re-deriving the reference with the corrected catalog as the one step not completed.
- **`mean_raw`** is the scoring method used throughout the paper. Switch to **`scanpy_corrected`** when comparing cells with very different transcriptome complexity (e.g., proliferating vs post-mitotic) — this prevents a baseline-shift artifact that can bias `mean_raw` scores.
- **Readiness vs activity:** module *readiness* (R, mean expression of the core machinery) and module *activity* (A, mean expression of TF target genes) are different things. They correlate at ρ ≈ 0.75 across modules, but CREB and NPAS4 are systematic exceptions because their activation is post-translational. The catalog flags `dissociation_risk: HIGH` for these — use both R and A when working with them.
- **The eigenspace is a fixed reference**; new datasets project into it without refitting. See [`docs/EIGENSPACE_EVOLUTION.md`](docs/EIGENSPACE_EVOLUTION.md) for the relationship to the first draft's 12-PC space. The axes are individually interpretable (PC1 inflammatory → terminal germline cells, PC2 quiescent → proliferative, PC3 signaling-dominated immune → metabolic epithelia), but convergence is a high-dimensional nearest-neighbour property, not a 2-D location: a raw PC1–PC2 projection recovers 8.2% core assignment for B-ALL against 60.8% in the full nine dimensions.
- **The capacity-floor predictor is upward-asymmetric.** A cell with `A_baseline > 4.5` for a module cannot ramp UP under stimulus. It CAN be actively suppressed DOWN by specific signaling perturbations (e.g., retinoic acid → UPR-ATF6). The predictor is a ceiling test, not a "this module won't change" claim.
- **What the criteria are.** The four operational criteria (a sensor, a transcription-factor output with targets, a negative-feedback element, receptor-field specificity) record *how the panel was curated*; they are not a verified property of the 206 genes as scored. A feedback element is named for 38 of 44 modules (six are marked unclear), and for ten of those the named element is not among the scored genes.

---

## Scientific background

The framework was developed across a series of papers (Spiro 2024-2026). This section briefly describes them so you understand what the tool implements and why.

- **Paper 3** — *Diverse cancers converge on one region of normal cell-state space in a cancer-blind coordinate system* ([Zenodo 10.5281/zenodo.20542130](https://doi.org/10.5281/zenodo.20542130)). Defines the 44-module catalog, the fixed nine-dimensional reference space built from 154 normal HPA cell types, and the convergence assay summarised above. It supersedes the first Paper 3 draft ("The perceptual eigenspace of the cell: 43 transcriptional modules organize cancer and aging"; 43 modules, 12-PC eigenspace, v0.1.0 of this tool): NPAS4 is the 44th module — a leave-one-out rebuild on 43 modules reproduces all nine components with eigenvector cosine ≥ 0.95 — and the final paper describes the axes differently and does not describe convergence as a single axis.
- **Paper 4** — *memory consolidation as neuronal perceptual cycle.* Adds NPAS4, a neuron-specific immediate-early signaling pathway. Results summarized at [mool32.github.io/memory-consolidation](https://mool32.github.io/memory-consolidation/).
- **Paper 4.1** — *perceptome on drugs.* Tested whether the framework can recover drug mechanisms from CMap LINCS L1000 transcriptional signatures. Closed with **6 validated mechanism-pathway recoveries** (MEKi/ERK, Proteasomei/UPR-PERK, CDKi/Cell Cycle, EGFRi/ERK, PI3Ki/PI3K-PTEN, IKKi/NF-κB) plus 3 stress-axis positive controls — and **11 pre-registered falsifications**, including drug-disease cosine matching in three independent formulations. The tool implements the validated operation (`pct.activity_layer_screen` + `pct.drug_anchors`) and explicitly does not support the falsified ones.
- **Paper 4.2** — *cancer transformation as a direction in capacity space* (legacy). Identified an 8-cell attractor cluster for 11–13 cancers and a per-module Δ vector (bundled as `attractor_v1.json`; P3 PASS in the Sun 2021 paired HCC cohort). Paper 3 extends and supersedes the cluster as the account of where cancers converge; the capacity-direction vector remains available.
- **Papers 4.3 - 4.8** — *substrate-series.* Six independent paradigms (immune cell activation, muscle hypertrophy, intestinal epithelium, organoid differentiation, organoid regeneration, organoid retinoid perturbation) testing the capacity-floor predictor. Closed with the predictor validated in active differentiation / acute remodeling contexts and a clear boundary: the architecture is NOT observed in steady-state homeostatic terminally-differentiated cells. Refined to upward-asymmetric (saturated cells cannot ramp UP, but specific perturbations can suppress DOWN).

Every operation has a documented validation status in [`docs/SCOPE.md`](docs/SCOPE.md), which records falsified and retracted claims with the same prominence as validated ones; the Paper 4.1 falsification audit trail is reproduced in the source tree (`killed_hypotheses.csv`).

## Tests

```bash
pip install pytest
pytest tests/ -v
# 205 passing in a few seconds
```

Coverage: catalog (and its agreement with the paper's Table S9 / S9b), scoring and gene aliases, standardized projection and calibrations, the convergence pipeline (assignment, quarantine, shape-only, within-patient test with the paper's exact statistical conventions, complexity matching, recurrence ranking), bundled paper tables (with the totals the paper itself states), statements of the paper that the bundled reference reproduces (Kaiser 9, 73.8% of variance, Fig. S1 stability tiers, NPAS4 leave-one-out), perceptivity (with biological regression checks against published a-priori predictions), attractor reference, validity scorecard, drug anchors + activity_layer_screen, end-to-end pipeline, reproducibility.

## Citation

If you use perceptome in your research, please cite the paper and the software.

**Paper:**

```bibtex
@misc{spiro2026convergence,
  author    = {Spiro, Theodor},
  title     = {Diverse cancers converge on one region of normal cell-state space in a cancer-blind coordinate system},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.20542130}
}
```

**Software (this version):**

```bibtex
@software{spiro2026perceptome,
  author    = {Spiro, Theodor},
  title     = {perceptome: Cellular Perception Analysis Toolkit},
  version   = {0.2.2},
  year      = {2026},
  publisher = {Zenodo},
  doi       = {10.5281/zenodo.20113710},
  url       = {https://github.com/mool32/perceptome}
}
```

**Concept DOI** (always resolves to the latest version): [10.5281/zenodo.20113468](https://doi.org/10.5281/zenodo.20113468)

For reproducibility, cite the **version DOI** of the specific release you used:
- v0.2.2: `10.5281/zenodo.20113710`
- v0.2.1: `10.5281/zenodo.20113516`

## Versioning

The latest release is **v0.2.2**. The `[Unreleased]` section of [`CHANGELOG.md`](CHANGELOG.md) holds the alignment with the final Paper 3: the standardized `pct.project` (which changes its output — `standardize=False` restores the old behaviour), reference calibrations, the `pct.convergence` subpackage and its bundled tables, gene-symbol handling, Table S9b criteria, and corrections of statements the paper retracts. See [`ROADMAP.md`](ROADMAP.md) for what comes next.

## Contributing

Contributions are welcome. Bug fixes and documentation are casual; new modules and reference vectors require the same discipline used to build the framework. See [`CONTRIBUTING.md`](CONTRIBUTING.md) for the full bar (in short: new canonical modules need the four operational criteria documented as in Table S9b — with `unclear` where a feedback element cannot be named — plus a published a-priori prediction reproduced as a regression test; falsified operations from Paper 4.1 will be rejected).

## License

MIT — see [`LICENSE`](LICENSE).

## Author

Theodor Spiro · ORCID [0009-0004-5382-9346](https://orcid.org/0009-0004-5382-9346) · tspiro@vaika.org · Vaika Inc., East Aurora, NY, USA

Issues and pull requests welcome on [GitHub](https://github.com/mool32/perceptome).
