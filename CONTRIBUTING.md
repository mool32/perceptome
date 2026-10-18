# Contributing to perceptome

Contributions are welcome. The bar varies by what you're contributing — bug fixes and documentation are casual, new scientific content (modules, reference vectors, drug anchors, validity controls) requires the same pre-registration discipline used to build the framework in the first place.

## Quick types

| Contribution | Bar | Process |
|---|---|---|
| **Bug fix** | Reproduce on `main`, write a test that fails before the fix and passes after | Standard PR |
| **Documentation** (README, QUICKSTART, SCOPE, docstrings) | Read for accuracy; don't introduce claims not backed by the code | Standard PR |
| **New utility function / convenience wrapper** | Doesn't change the framework's scientific meaning; covered by ≥1 unit test | Standard PR |
| **Performance optimization** | Doesn't change numerical output; benchmark in PR description; covered by existing tests | Standard PR |
| **New scoring method** | Document scope (when to use vs default), benchmark against `mean_raw` on the bundled tutorial data | PR with discussion in an issue first |
| **New module added to the catalog** | See "New canonical module" below — high bar | PR + issue + reference |
| **New reference vector** (disease, aging, attractor variant, drug anchor) | See "New reference vector" below | PR + issue + script |
| **A statement about the convergence result** (docs, docstrings, notebooks) | Must match what Paper 3 supports; see "Statements about convergence" below | Standard PR |

## Standard PR process

1. Open an issue describing the change before opening the PR (especially for anything beyond a bug fix). Saves both of us time if the direction isn't going to be accepted.
2. Fork the repo, branch off `main`, make your changes.
3. Run the test suite and make sure it passes (it uses only the bundled data — no network, no downloads):
   ```bash
   pytest tests/ -v
   ```
4. If you added functionality, add a test. If you changed numerical behavior, update the relevant biological regression checks in `tests/test_integration.py` and explain in the PR why the new value is correct.
5. Run linting / formatting (project uses no formatter currently — match the existing style).
6. Update `CHANGELOG.md` under an `[Unreleased]` section.
7. Open the PR with a clear description: what, why, and what you tested.

## High bar: new canonical module

Adding a 45th (or 46th, ...) module to the catalog is a scientific change, not a coding change. In Paper 3 a module is "a self-contained regulatory unit specified by four operational criteria", and a candidate is kept only if it represents a biologically distinct regulatory input or response, if its activity can be inferred from a defined transcriptional output rather than from the expression of a single pathway component, and if its sensor, transcriptional regulator and target-gene response are characterized well enough to fix the direction of pathway activity. The four criteria:

1. **Signal-sensing component** — an identifiable upstream input (ligand-receptor, metabolite, mechanical, …)
2. **Transcription-factor output** — a transcription factor (or family) with documented target genes, measurable in scRNA-seq
3. **Cyclical negative-feedback element** — a feedback element induced by the module's own activity
4. **Receptor-field specificity** — the class of input the sensor reads, written down. A candidate that substantially duplicates a module already in the catalog, lacks a defined transcriptional output, or acts mainly post-transcriptionally is dropped.

**Document it the way Table S9b of Paper 3 does** (`pct.module_criteria()` shows the format): the sensor, the TF output, the feedback element and the receptor field, in prose, with a confidence grade. Three rules come from how the paper audits its own table:

- **Say "unclear" rather than filling a criterion in.** Six of the 44 modules carry an unclear feedback element.
- **Check the named feedback element against the genes you actually score.** For ten modules it is not among them, so the score cannot read the feedback arm; record that.
- **Never write a criterion back retrospectively without saying so.** The receptor field is blank for the 16 modules annotated at catalog construction, which the paper reports as a gap in the record rather than filling it in.

Be clear in the PR about what the criteria are: the paper's own audit treats them as a record of how the panel was curated, not as a verified property of the scored genes (three curators given the criteria alone built materially different catalogs).

Plus, infrastructure for the new module:

- **Gene panel** — `core_genes` (3-5 sensor cascade + TF), `activity_genes` (5-15 TF target genes), optional `sensor_genes` / `cascade_genes` / `tf_genes` / `feedback_genes`
- **Gene symbols** — current HGNC symbols, with the coverage checked on the HPA reference and on at least one cohort (`pct.score_modules(...)["coverage"]` and `["missing_genes"]`). The paper lost three catalog genes to renamed symbols and to a loader that stripped `-<digits>` suffixes (Note S-Q); record any renamed symbol in `perceptome/catalog/aliases.py`.
- **Recompute scripts updated** — `scripts/01_build_catalog_v03.py`, `scripts/02_build_hpa_perceptivity.py` (HPA reference now N × 45), `scripts/03_build_eigenspace.py` (eigenspace rebuilt from N × 45 R matrix). Then `scripts/04_build_attractor_reference.py`, `scripts/06_recompute_disease_vectors.py` and `scripts/07_recompute_aging_reference.py` rerun on the updated eigenspace, and `scripts/13_build_single_cell_calibration.py` (the bundled single-cell calibration, accepted only if it still reproduces Table S7).
- **Does the space survive the change?** — rebuild the eigenspace and report the Kaiser count and the per-PC eigenvector cosine against the previous reference, as the paper's NPAS4 leave-one-out does (Note S2); `scripts/15_paper_reproduction_report.py` shows which published numbers the new reference still reproduces. The paper's tables bundled in `perceptome/convergence/data/` belong to the 44-module catalog; say so in `docs/PAPER_MAP.md` rather than leaving them to look current.
- **At least one published a-priori prediction reproduced** — add a regression test to `tests/test_integration.py` of the form "starting cell C, module M, expected `capacity_floor` = X, observed = X" with citation to the empirical study that pre-registered the prediction.
- **Updated `docs/SCOPE.md`** — record the validation evidence + scope conditions for the new module
- **CHANGELOG entry** — under "Added" with explicit version bump rationale (catalog change = at least minor bump, e.g. v0.2 → v0.3)

If you don't have a published a-priori prediction yet — that's fine, the work is at the empirical-research stage, not the tool-contribution stage. The catalog stays at 44 until the validation exists.

## High bar: new reference vector

Adding a disease, aging, attractor, or drug reference vector means another column in the `compare_to_references` cosine output. The bar is similar to a module:

- **Source documented** — what scRNA-seq dataset, what comparison (disease vs control, old vs young, tumor vs normal), what cell-group definitions
- **Reproducible build script** — added to `scripts/`, follows the pattern of `06_recompute_disease_vectors.py` or `07_recompute_aging_reference.py`
- **Standardized coordinates** — store the per-module delta and derive the eigenspace vector with `pct.project_shift` (module delta divided by the reference SD, then the fixed loadings), and set `eigenspace_convention: standardized` in the file. Do not project a difference of raw scores through the loadings: that is the pre-0.3 convention and its cosines are not comparable with the bundled vectors.
- **At least one validation test** — vector loads, has expected sign on a key module
- **Avoid double-counting** — don't add a vector that's near-duplicate of an existing one (cosine > 0.85 with an existing reference); explain how it's complementary

## High bar: statements about convergence

Claims about the convergence result are held to what Paper 3 supports. Read [`docs/SCOPE.md`](docs/SCOPE.md) first: the region is the robust finding and the anchor is catalog-conditional; the within-patient increase is the unit of evidence (the patient, not the cell); convergence is a static endpoint, not a measured trajectory; and a list of statements the paper retracted (baseline gating, fold-change framing, the cGAS-STING reading, the polyploidy argument) must not come back in code, docstrings or documentation. A change to `pct.convergence` that alters a number the paper prints needs the paper's value, the new value and the reason in the PR; `tests/test_paper_reference.py` and `tests/test_convergence.py` hold those numbers.

## What gets rejected

The framework has been built through a pre-registered research program with 11 falsifications recorded in `killed_hypotheses.csv` and `PAPER4_1_CANONICAL_RESULTS.md`. The tool deliberately does not implement operations that have been tested and shown to fail. Examples of contributions that will be rejected:

- **A `drug_disease_cosine_match()` function**, even with caveats. Three independent formulations were tested and falsified in Paper 4.1; reintroducing it in the tool would invite users to do exactly the analysis the data refused.
- **A `tf_autoregulation_layer()` predictor**. Tested in Paper 4.1 v3 (8 drug classes, K=1/8 PASS, binomial P(K≥1)=0.57 = pure noise). Killed.
- **A `clean_signature_score()` linked to FDA approval**. Tested in Paper 4.1 program-engagement v1 (ΔR² = +0.002, interaction wrong direction, selectivity captures all signal). Killed.
- **Any `predict_drug_efficacy()` claim** without per-cancer per-cell-line replication. The 6 surviving Paper 4.1 anchors are mechanism-pathway recoveries on CMap, not efficacy predictors.

If you want to revisit one of the 11 falsified hypotheses with new data or a new design, the path is: write a fresh pre-registration document with locked thresholds and scramble nulls, run it, get a PASS, document in a published paper, then propose the tool addition based on the published evidence. Not the other way around.

## Code style

No automated formatter currently. Match the existing style:

- Module docstrings explain *what the module is for* and *what's reserved for later versions*
- Function docstrings: NumPy style, short Parameters / Returns sections, no exhaustive type narratives
- Comments only when the *why* is non-obvious (a workaround, a paper-specific constant, a calibration choice). Don't comment what the code obviously does.
- No emoji in source files (README and other docs OK if the user-facing tone calls for it)
- Tests should be small, named after what they verify, and prefer biological assertions over implementation assertions

## Testing

```bash
# Full suite (a few seconds on a laptop; no network needed)
pytest tests/ -v

# Single file while iterating
pytest tests/test_perceptivity.py -v

# Coverage
pytest tests/ --cov=perceptome --cov-report=term-missing
```

The test suite includes biological regression checks (paper4.4 / 4.5 / 4.7 a-priori predictions, attractor direction signs, HPA reference v1.4 spot-checks) and, since the Paper 3 alignment, the numbers of that paper that the bundled data can reproduce (`tests/test_paper_reference.py`: component count, variance, stability tiers, gene lists, Table S9b counts, the totals of the transcribed tables) and the convergence workflow on cells synthesized from the bundled reference (`tests/test_convergence.py`, `tests/test_pipeline_qc_aliases.py`). If a contribution changes one of these numbers, the PR must justify why and reference the empirical work that supports the new value.

## Getting in touch

- Bug reports + feature requests: [GitHub Issues](https://github.com/mool32/perceptome/issues)
- Scientific questions about the framework: open a Discussion or contact tspiro@vaika.org
- Security issues: tspiro@vaika.org (do not file public issue)

By contributing, you agree your contributions are licensed under the same MIT license as the rest of the project.
