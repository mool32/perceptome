# perceptome

[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](LICENSE)
[![Python: ≥3.8](https://img.shields.io/badge/Python-≥3.8-blue.svg)](https://www.python.org)
[![Tests: 73/73](https://img.shields.io/badge/tests-73%2F73%20passing-green.svg)](#tests)
<!-- ZENODO_DOI_BADGE_PLACEHOLDER -->
[![DOI: pending](https://img.shields.io/badge/DOI-pending-lightgrey.svg)](https://github.com/mool32/perceptome/releases/tag/v0.2.1)

A Python toolkit for analyzing single-cell RNA-seq data through the lens of **44 cellular signaling pathways** (NF-κB, mTOR, UPR, p53, nuclear receptors, …) treated as a **perceptual system** — the machinery a cell uses to sense and respond to its environment.

For each cell or cell type, perceptome answers three questions other tools don't directly answer:

1. **Where does this cell sit in signaling-module space?** A 9-dimensional coordinate system derived from 154 normal human cell types (Human Protein Atlas), so any cell — yours, a tumor cell, a drug-treated cell — can be located on the same map.
2. **What can this cell do under perturbation?** A capacity-floor predictor: given a cell type and a pathway, will it ramp up under stimulus or is it already saturated? Validated across multiple experimental paradigms (intestinal homeostasis, organoid regeneration, muscle hypertrophy).
3. **Did the measurement reach biology, or hit a confound?** A built-in scorecard with three null controls (random-gene panel, housekeeping, cell cycle) that catches the common artifacts in perturbation analyses before they corrupt downstream interpretation.

The framework is the result of a multi-year, pre-registered research program. The tool implements **only what those studies validated** and explicitly excludes operations that were tested and falsified — see [`docs/SCOPE.md`](docs/SCOPE.md).

---

## Try it in 5 minutes

- 📓 [`examples/tutorial.ipynb`](examples/tutorial.ipynb) — full executable walkthrough on PBMC3K (the standard 10x Genomics demo dataset). Auto-downloads, runs in <2 minutes, ships with output cells populated. **Best starting point.**
- 📖 [`docs/QUICKSTART.md`](docs/QUICKSTART.md) — same walkthrough in markdown, with extra notes on common pitfalls.

```bash
pip install git+https://github.com/mool32/perceptome.git
```

```python
import scanpy as sc
import perceptome as pct

adata = sc.read_h5ad("your_data.h5ad")              # log-normalized
scores = pct.score_modules(adata)["scores"]         # 44 modules per cell
coords = pct.project(scores)["coordinates"]         # 9-PC eigenspace
refs = pct.compare_to_references(coords)            # cosines vs cancer / disease / aging
```

That's the geometry workflow. There are three more — see the tutorial.

---

## What problem does this solve?

Most single-cell analyses operate at the gene level (DEGs, marker genes, GSEA on gene sets) or at the cluster level (cell type labeling, trajectory inference). Both are useful. Neither directly answers questions like:

- *"My drug shifts cells in the dataset — does that shift align with the direction tumors take during transformation?"*
- *"My cells of type X — can they ramp pathway Y under stimulus, or is that pathway already at ceiling?"*
- *"My perturbation effect on module Z — is it real biology or a baseline artifact from the cells being more proliferative?"*

perceptome works one level up: **44 pathways, treated as the cell's signaling repertoire**, projected into a fixed coordinate system where any cell can be compared to any reference — and a capacity layer that predicts the ramp room of each pathway in each cell type before you run the experiment.

## When to use it / when not

**Good fit:**
- You have scRNA-seq (or bulk RNA-seq) and want a low-dimensional, interpretable view of pathway activity per cell type
- You're designing or interpreting a perturbation experiment (drug, knockdown, stimulus) and want a-priori predictions of which pathways have ramp room in which cell types
- You're working on cancer transformation, disease perturbation, or aging — and want to test alignment with directions derived from independent reference cohorts
- You need a built-in validity check before trusting perturbation effects

**Not a good fit:**
- Single-gene differential expression — use `scanpy.tl.rank_genes_groups` or DESeq2 directly
- Cell type identification or label transfer — use scanpy / CellTypist / Symphony
- Trajectory inference — use Slingshot / scFates / palantir
- Drug-disease repositioning via geometric similarity — this was tested in three independent ways and falsified; the tool deliberately doesn't support it (see [`docs/SCOPE.md`](docs/SCOPE.md))
- Phosphorylation cascade readout — the framework operates on transcription, not protein phosphorylation
- Clinical decision support — research tool only, no clinical validation

---

## What's in the box

| Subpackage | Function | What it answers |
|---|---|---|
| **catalog** | `pct.list_modules()`, `pct.get_genes(module, gene_set)` | What modules and gene panels are defined? |
| **score** | `pct.score_modules`, `pct.score_readiness`, `pct.score_activity` | Module activity per cell from expression |
| **eigenspace** | `pct.project`, `pct.rebuild` | Project into the 9-PC reference perceptual space |
| **perceptivity** | `pct.compute_perceptivity`, `pct.predict_engagement`, `pct.capacity_floor` | Capacity, headroom, saturation per cell type × module |
| **compare** | `pct.compare_conditions`, `pct.compare_to_references`, `pct.infrastructure_regime`, `pct.divergence_score` | Disease vs control deltas, regime classification |
| **reference** | `pct.load_attractor_direction`, `pct.attractor_alignment` | Cancer transformation direction in capacity space |
| **validity** | `pct.validate_perturbation`, `pct.scorecard` | Random-200 / housekeeping / cell-cycle nulls |
| **drugs** | `pct.drug_anchors`, `pct.activity_layer_screen` | 9 validated mechanism-pathway anchors (narrow scope) |
| **network** | `pct.compute_network`, `pct.module_heterogeneity` | Module co-variation + rich-club analysis |

## Bundled reference data

- **44-module catalog** with core / activity / sensor / cascade / TF / feedback gene lists per module
- **154 × 44 cell-type-by-module reference matrix** computed from the Human Protein Atlas, with R (readiness), A (activity), C (capacity), and headroom precomputed
- **9-PC reference eigenspace** built from the same 154 × 44 matrix; PC1-PC6 stable per 100 bootstraps
- **Cancer transformation direction** — 8-cell-type attractor cluster + locked Δ vector in module + eigenspace coordinates
- **Disease reference vectors** for rheumatoid arthritis, Alzheimer's, idiopathic pulmonary fibrosis, and diabetic kidney disease (all recomputed from raw scRNA-seq on the v0.2 eigenspace)
- **Aging reference vectors** — inflammaging direction (blood old−young) and collapse direction (bone marrow old−young)
- **9-anchor drug reference table** with full validation chain
- **Validity panels**: housekeeping (Eisenberg & Levanon 2013), cell cycle (Tirosh 2016)

## Methodology notes

- **`mean_raw`** is the default scoring method; safe across datasets. Switch to **`scanpy_corrected`** when comparing cells with very different transcriptome complexity (e.g., proliferating vs post-mitotic) — this prevents a baseline-shift artifact that can bias `mean_raw` scores.
- **Readiness vs activity:** module *readiness* (R, mean expression of the core machinery) and module *activity* (A, mean expression of TF target genes) are different things. They correlate at ρ ≈ 0.75 across modules, but CREB and NPAS4 are systematic exceptions because their activation is post-translational. The catalog flags `dissociation_risk: HIGH` for these — use both R and A when working with them.
- **Eigenspace is a fixed reference**; new datasets project into it without refitting. See [`docs/EIGENSPACE_EVOLUTION.md`](docs/EIGENSPACE_EVOLUTION.md) for the relationship between the v0.1 (12-PC) and v0.2 (9-PC) reference spaces if you're cross-referencing earlier results.
- **The capacity-floor predictor is upward-asymmetric.** A cell with `A_baseline > 4.5` for a module cannot ramp UP under stimulus. It CAN be actively suppressed DOWN by specific signaling perturbations (e.g., retinoic acid → UPR-ATF6). The predictor is a ceiling test, not a "this module won't change" claim.

---

## Scientific background

The framework was developed across a series of papers (Spiro 2024-2026, in preparation). They are not yet preprints; this section briefly describes them so you understand what the tool implements and why.

- **Paper 3** — *the perceptome eigenspace.* Introduces the 43-module catalog and the 12-PC eigenspace derived from 154 Human Protein Atlas cell types. Documents perception breadth as PC1, the cancer convergence axis on PC4, and the four-quadrant infrastructure regime classifier (supply chain / firefighting / collapse / unsupported). v0.1 of this tool accompanied Paper 3.
- **Paper 4** — *memory consolidation as neuronal perceptual cycle.* Adds NPAS4 as the 44th canonical module (a neuron-specific immediate-early signaling pathway). v0.2 of this tool incorporates this addition; the eigenspace was rebuilt accordingly (now 9 PCs from 154 × 44).
- **Paper 4.1** — *perceptome on drugs.* Tested whether the framework can recover drug mechanisms from CMap LINCS L1000 transcriptional signatures. Closed with **6 validated mechanism-pathway recoveries** (MEKi/ERK, Proteasomei/UPR-PERK, CDKi/Cell Cycle, EGFRi/ERK, PI3Ki/PI3K-PTEN, IKKi/NF-κB) plus 3 stress-axis positive controls — and **11 pre-registered falsifications**, including drug-disease cosine matching in three independent formulations. The tool implements the validated operation (`pct.activity_layer_screen` + `pct.drug_anchors`) and explicitly does not support the falsified ones.
- **Paper 4.2** — *cancer transformation as direction in capacity space.* Identified an 8-cell-type "attractor" cluster of normal cell states that 11 cancer types from 11 organ systems converge toward (mean alignment 42% vs 5% chance baseline). Independently replicated on Sun et al. 2021 paired hepatocellular carcinoma (cosine alignment > +0.20 in 4/6 cell types within paired patients). The bundled `attractor_v1.json` contains this direction.
- **Papers 4.3 - 4.8** — *substrate-series.* Six independent paradigms (immune cell activation, muscle hypertrophy, intestinal epithelium, organoid differentiation, organoid regeneration, organoid retinoid perturbation) testing the capacity-floor predictor. Closed with the predictor validated in active differentiation / acute remodeling contexts and a clear boundary: the architecture is NOT observed in steady-state homeostatic terminally-differentiated cells. Refined to upward-asymmetric (saturated cells cannot ramp UP, but specific perturbations can suppress DOWN).

The pre-registration discipline is part of the tool itself: every operation has a documented validation status in [`docs/SCOPE.md`](docs/SCOPE.md), and the Paper 4.1 falsification audit trail is reproduced in the source tree (`killed_hypotheses.csv`).

When the preprints are posted, this section will link to them directly.

## Tests

```bash
pip install pytest
pytest tests/ -v
# 73/73 passing in <2 seconds
```

Coverage: catalog, scoring, perceptivity (with biological regression checks against published a-priori predictions), eigenspace projection, attractor reference, validity scorecard, drug anchors + activity_layer_screen, end-to-end pipeline, reproducibility.

## Citation

A Zenodo DOI is being assigned for v0.2.0. Once available, citation will be:

```bibtex
@software{spiro2026perceptome,
  author = {Spiro, Theodor},
  title = {perceptome: Cellular Perception Analysis Toolkit},
  version = {0.2.0},
  year = {2026},
  doi = {10.5281/zenodo.XXXXXXX},
  url = {https://github.com/mool32/perceptome}
}
```

The DOI badge will be added here once Zenodo issues it.

## Versioning

This is **v0.2.0** — first widely-published release. Major additions vs v0.1.0:

- 44-module catalog (NPAS4 added)
- 9-PC eigenspace recomputed
- New `pct.perceptivity` subpackage (capacity layer with two-factor framework)
- New `pct.validity` subpackage (perturbation-analysis scorecard)
- New `pct.reference.attractor` (cancer transformation direction)
- New `pct.drugs` subpackage (narrow validated drug screen)
- Disease + aging reference vectors recomputed on the new eigenspace
- Tutorial notebook + 3 documentation files

See [`CHANGELOG.md`](CHANGELOG.md) for the full diff and [`ROADMAP.md`](ROADMAP.md) for v0.3+ direction.

## License

MIT — see [`LICENSE`](LICENSE).

## Author

Theodor Spiro · ORCID [0009-0004-5382-9346](https://orcid.org/0009-0004-5382-9346) · tspiro@vaika.org · Vaika Inc., East Aurora, NY, USA

Issues and pull requests welcome on [GitHub](https://github.com/mool32/perceptome).
