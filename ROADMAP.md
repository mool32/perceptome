# perceptome Roadmap

## v0.2.0 — Capacity layer + restructured architecture (2026-05) [CURRENT]

Released alongside the substrate-series papers (4.3-4.8) and Paper 4.2 capacity extension. Adds the perceptivity / capacity layer, validity scorecard, attractor reference, focused 9-anchor drug operation. See [`CHANGELOG.md`](CHANGELOG.md) for full diff.

## v0.2.1 — Validation expansion [NEXT]

The v0.2.0 release closed the major v0.1 → v0.3 transition (catalog, eigenspace, perceptivity layer, attractor, validity, drugs, disease/aging vectors all rebuilt on the new basis). v0.2.1 focus shifts to validation breadth and quality-of-life:

- **Tabula Sapiens cross-atlas validation** for the perceptivity metric (Criterion 2 in metric design v1.4 — deferred from v0.2.0).
- **Statistical sanity tests** for `activity_layer_screen` (verify p-distribution under known null).
- **Edge-case coverage** in tests (empty AnnData, single cell, single module, all-NaN expression).
- **Cross-version regression fixture** (known-good output on a small open dataset for drift detection).
- **Disease vector expansion**: add Type 2 diabetes, Parkinson's, ALS as additional reference points if raw data accessible.

## v0.3 — Tissue-specific module extensions [PLANNED]

Per Paper 4 ROADMAP: tissue-specific modules that pass P1–P4 criteria but are not pan-cellular (currently flagged via `pan_cellular: bool` and `tissue_bias: list[str]` in the catalog).

- HNF4α (liver-specific)
- MYOD (muscle-specific)
- Tissue-specific catalog entries with explicit pan-cellular vs tissue-specific flags
- Conditional eigenspace projection: project against pan-cellular eigenspace OR tissue-specific extended eigenspace (per cell-type-class context)

## v0.4 — Two-factor framework Factor 2 [PLANNED]

Currently `predict_engagement()` returns Factor 1 (capacity) only and `predicted_magnitude='unknown'` always. Factor 2 (operation intensity) requires:

- Operation taxonomy (acute differentiation, hypertrophy, regeneration, terminal differentiation, ...)
- Per-operation intensity calibration on substrate-series data (paper4.4-4.8 + paper4)
- Factor 2 estimator from observed substrate-series effect sizes
- Joint magnitude prediction: `predicted_magnitude = capacity × operation_intensity`

This closes the metric design v1.4 framework and lets the tool answer "how big a ramp" instead of just "ramp possible / blocked".

## v0.5 — PyPI publication + stable API [VISION]

- API guarantees: scoring, projection, perceptivity, regime classification stable across patch versions
- Publish to PyPI (`pip install perceptome`)
- Pre-computed canonical projections for major public datasets (HCA, Tabula Sapiens, GTEx, CELLxGENE Census)
- Cross-species support: mouse ortholog mapping handled internally (current ~95% concordance via uppercase mapping)

## v1.0 — Dynamics on the eigenspace [VISION]

Long-term: dynamical models on the perceptome eigenspace potential landscape (trajectories instead of static positions). Requires the substrate-series Factor 2 closure and capacity-trajectory paired-cohort data (e.g., Sun 2021-style designs across multiple cancers and disease contexts).

---

## Companion paper schedule (independent of tool versions)

- **Paper 3** — perceptome eigenspace + cancer + aging — accompanied v0.1.0 (43 modules, 12 PCs)
- **Paper 4** — Memory consolidation as neuronal perceptual cycle (NPAS4 as 44th module) — drove v0.2.0 catalog expansion
- **Paper 4.1** — Perceptome on drugs: 6 validated rescues + 11 falsifications + 5 SC scope conditions (closed 2026-05-09) — drove the v0.2.0 drug subpackage rewrite
- **Paper 4.2** — Convergent perceptual identity in cancers across organ systems (capacity-direction P3 PASS in Sun 2021 paired HCC) — drove `pct.reference.attractor`
- **Papers 4.3-4.8** — Substrate-series (immune, muscle, epithelial, organoid, organoid_regen, organoid_holdout) — closed the capacity-floor predictor (upward-asymmetric); drove `pct.perceptivity.floor`

## How to influence the roadmap

- File a GitHub issue for feature requests, bugs, or scientific questions about the framework.
- Pull requests welcome. For larger contributions, please open an issue first to coordinate.
- Citation graph is the strongest signal: papers using perceptome v0.2 to validate predictions in new biological systems directly inform v0.3+ priorities.

---

*Maintained by Theodor Spiro (Vaika Inc.). See [README](README.md) and [LICENSE](LICENSE) for details.*
