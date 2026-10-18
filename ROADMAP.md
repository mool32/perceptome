# perceptome Roadmap

## Released

- **v0.2.2** (2026-05-10) — discoverability and community patch: `CONTRIBUTING.md`, hero figure, two recipes.
- **v0.2.1** (2026-05-10) — documentation accessibility patch and the executable tutorial notebook.
- **v0.2.0** (2026-05-10) — capacity layer and restructured architecture: perceptivity / capacity metric, validity scorecard, attractor reference, focused 9-anchor drug operation, 44-module catalog and the 9-PC reference. See [`CHANGELOG.md`](CHANGELOG.md).

## v0.3 — Paper 3 alignment [NEXT — implemented under `[Unreleased]` in the changelog]

Aligns the tool with the final Paper 3, Spiro (2026), *Diverse cancers converge on one region of normal cell-state space in a cancer-blind coordinate system* ([Zenodo 10.5281/zenodo.20542130](https://doi.org/10.5281/zenodo.20542130)), and corrects the statements it retracts. In place: the standardized projection into the nine-dimensional reference space, reference calibrations, `pct.convergence` (assignment to 152 normal cell types, region occupancy, within-patient test), the paper's Supporting-Information tables, gene-symbol handling, the operational criteria of Table S9b, and the documentation mapping every claim to the code ([`docs/PAPER_MAP.md`](docs/PAPER_MAP.md)).

Before it is tagged:

- **Maintainer review** of [`docs/PAPER_MAP.md`](docs/PAPER_MAP.md): what the bundled data reproduce and what they do not (the cGAS-STING row of Table S7; the 0.99 among the origins), and the observations for the author at its end. The single-cell-scaled calibration is bundled and verified against Table S7; confirm that shipping the 154 × 44 summary of the HPA table (CC BY-SA 3.0) is acceptable.
- **Choose the version number, tag, and mint the Zenodo version DOI**, then update `CITATION.cff`. The default output of `pct.project` changes, which under 0.x semantics argues for a minor bump. The package still reports 0.2.2.

## Follow-ups the paper leaves open [PLANNED]

- **Re-derive the reference with the corrected catalog.** Three genes were lost to renamed symbols or a loader (Note S-Q); the paper states the corrected catalog "has not been re-derived end to end". Then switch `resolve_aliases` on by default and rebuild the eigenspace, calibrations and reference vectors.
- **A copy-number-based malignant floor** for the convergence workflow (the paper defers it "to production"): malignant-cell purity is the documented failure mode of the within-patient test.
- **The rule that selects the reported modules in the occupant comparison** (Fig. 6; `pct.occupant_contrast` returns P and a Benjamini-Hochberg q and leaves the threshold open), and replication outside colorectal.
- **A bundled normal single-cell reference for the matched-averaging null** (Tabula Sapiens or the Human Cell Atlas, which also serves the cross-atlas validation below) with a definition of a coherent micro-population, and the specification curve of the shape contraction once its variants are defined (`pct.contraction_null` is in place in generic form).
- **Cohorts the paper names as missing:** esophageal dysplasia (the pre-invasive prediction failed in non-dysplastic Barrett's metaplasia), distant rather than nodal metastases, a normal breast arm with more than one sample.
- **Recompute the empirical convergence direction** from the 25-cancer cohort; the bundled attractor direction is still the Paper 4.2 estimate (11–13 cancers).
- **A larger example gallery** (further cohorts through `pct.convergence`).

## v0.4 — Validation expansion [PLANNED] (formerly planned as v0.2.1)

- **Tabula Sapiens cross-atlas validation** for the perceptivity metric (Criterion 2 in metric design v1.4 — deferred from v0.2.0).
- **Statistical sanity tests** for `activity_layer_screen` (verify p-distribution under known null).
- **Edge-case coverage** in tests (empty AnnData, single cell, single module, all-NaN expression).
- **Cross-version regression fixture** (known-good output on a small open dataset for drift detection).
- **Disease vector expansion**: add Type 2 diabetes, Parkinson's, ALS as additional reference points if raw data accessible.

## v0.5 — Tissue-specific module extensions [PLANNED] (formerly v0.3)

Per Paper 4 ROADMAP: tissue-specific modules that pass P1–P4 criteria but are not pan-cellular (currently flagged via `pan_cellular: bool` and `tissue_bias: list[str]` in the catalog).

- HNF4α (liver-specific)
- MYOD (muscle-specific)
- Tissue-specific catalog entries with explicit pan-cellular vs tissue-specific flags
- Conditional eigenspace projection: project against pan-cellular eigenspace OR tissue-specific extended eigenspace (per cell-type-class context)

New modules are documented against the four operational criteria in the form of Table S9b of Paper 3 (`pct.module_criteria()`); see [`CONTRIBUTING.md`](CONTRIBUTING.md).

## v0.6 — Two-factor framework Factor 2 [PLANNED] (formerly v0.4)

Currently `predict_engagement()` returns Factor 1 (capacity) only and `predicted_magnitude='unknown'` always. Factor 2 (operation intensity) requires:

- Operation taxonomy (acute differentiation, hypertrophy, regeneration, terminal differentiation, ...)
- Per-operation intensity calibration on substrate-series data (paper4.4-4.8 + paper4)
- Factor 2 estimator from observed substrate-series effect sizes
- Joint magnitude prediction: `predicted_magnitude = capacity × operation_intensity`

This closes the metric design v1.4 framework and lets the tool answer "how big a ramp" instead of just "ramp possible / blocked".

## v0.7 — PyPI publication + stable API [VISION] (formerly v0.5)

- API guarantees: scoring, projection, perceptivity, regime classification stable across patch versions
- Publish to PyPI (`pip install perceptome`)
- Pre-computed canonical projections for major public datasets (HCA, Tabula Sapiens, GTEx, CELLxGENE Census)
- Cross-species support: mouse ortholog mapping handled internally (current ~95% concordance via uppercase mapping)

## v1.0 — Dynamics on the eigenspace [VISION]

Long-term: dynamical models on the perceptome eigenspace potential landscape (trajectories instead of static positions). Requires the substrate-series Factor 2 closure and capacity-trajectory paired-cohort data (e.g., Sun 2021-style designs across multiple cancers and disease contexts). Paper 3 states that its data are snapshots and that a static convergent endpoint is not a demonstrated attractor; time-resolved paired data are what would change that.

---

## Companion paper schedule (independent of tool versions)

- **Paper 3** — the first draft ("The perceptual eigenspace of the cell", 43 modules, 12 PCs) accompanied v0.1.0. The final version (Spiro 2026, [10.5281/zenodo.20542130](https://doi.org/10.5281/zenodo.20542130): 44 modules, nine-dimensional reference space, the convergence assay) is what v0.3 implements.
- **Paper 4** — Memory consolidation as neuronal perceptual cycle (NPAS4 as 44th module) — drove v0.2.0 catalog expansion
- **Paper 4.1** — Perceptome on drugs: 6 validated rescues + 11 falsifications + 5 SC scope conditions (closed 2026-05-09) — drove the v0.2.0 drug subpackage rewrite
- **Paper 4.2** — Convergent perceptual identity in cancers across organ systems (capacity-direction P3 PASS in Sun 2021 paired HCC) — drove `pct.reference.attractor`; its eight-cell cluster is the earlier estimate of the convergence target that Paper 3 (the two-cell core and the region around it) supersedes
- **Papers 4.3-4.8** — Substrate-series (immune, muscle, epithelial, organoid, organoid_regen, organoid_holdout) — closed the capacity-floor predictor (upward-asymmetric); drove `pct.perceptivity.floor`

## How to influence the roadmap

- File a GitHub issue for feature requests, bugs, or scientific questions about the framework.
- Pull requests welcome. For larger contributions, please open an issue first to coordinate.
- Citation graph is the strongest signal: papers using perceptome to validate predictions in new biological systems directly inform the priorities above.

---

*Maintained by Theodor Spiro (Vaika Inc.). See [README](README.md) and [LICENSE](LICENSE) for details.*
