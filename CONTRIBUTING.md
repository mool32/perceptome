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

## Standard PR process

1. Open an issue describing the change before opening the PR (especially for anything beyond a bug fix). Saves both of us time if the direction isn't going to be accepted.
2. Fork the repo, branch off `main`, make your changes.
3. Run the test suite and ensure all 73 tests still pass:
   ```bash
   pytest tests/ -v
   ```
4. If you added functionality, add a test. If you changed numerical behavior, update the relevant biological regression checks in `tests/test_integration.py` and explain in the PR why the new value is correct.
5. Run linting / formatting (project uses no formatter currently — match the existing style).
6. Update `CHANGELOG.md` under an `[Unreleased]` section.
7. Open the PR with a clear description: what, why, and what you tested.

## High bar: new canonical module

Adding a 45th (or 46th, ...) module to the catalog is a scientific change, not a coding change. The four perceptual criteria from the framework's foundation must hold:

1. **Sensor present** — the module is activated by an identifiable upstream signal (ligand-receptor, metabolite, mechanical, …)
2. **TF output** — the module has a primary transcription factor (or factor family) whose targets are measurable in scRNA-seq
3. **Cyclical feedback** — the module has a negative feedback loop closing within minutes-to-hours of activation
4. **Receptor field specificity** — the module is not a downstream consequence of one of the existing 44 modules

Plus, infrastructure for the new module:

- **Gene panel** — `core_genes` (3-5 sensor cascade + TF), `activity_genes` (5-15 TF target genes), optional `sensor_genes` / `cascade_genes` / `tf_genes` / `feedback_genes`
- **Recompute scripts updated** — `scripts/01_build_catalog.py`, `scripts/02_build_hpa_perceptivity.py` (HPA reference now N × 45), `scripts/03_build_eigenspace.py` (eigenspace rebuilt from N × 45 R matrix). The `scripts/04_build_attractor_reference.py` and `scripts/06_recompute_disease_vectors.py` will need to rerun on the updated eigenspace.
- **At least one published a-priori prediction reproduced** — add a regression test to `tests/test_integration.py` of the form "starting cell C, module M, expected `capacity_floor` = X, observed = X" with citation to the empirical study that pre-registered the prediction.
- **Updated `docs/SCOPE.md`** — record the validation evidence + scope conditions for the new module
- **CHANGELOG entry** — under "Added" with explicit version bump rationale (catalog change = at least minor bump, e.g. v0.2 → v0.3)

If you don't have a published a-priori prediction yet — that's fine, the work is at the empirical-research stage, not the tool-contribution stage. The catalog stays at 44 until the validation exists.

## High bar: new reference vector

Adding a disease, aging, attractor, or drug reference vector means another column in the `compare_to_references` cosine output. The bar is similar to a module:

- **Source documented** — what scRNA-seq dataset, what comparison (disease vs control, old vs young, tumor vs normal), what cell-group definitions
- **Reproducible build script** — added to `scripts/`, follows the pattern of `06_recompute_disease_vectors.py` or `07_recompute_aging_reference.py`
- **At least one validation test** — vector loads, has expected sign on a key module
- **Avoid double-counting** — don't add a vector that's near-duplicate of an existing one (cosine > 0.85 with an existing reference); explain how it's complementary

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
# Full suite (1.9 sec on a laptop)
pytest tests/ -v

# Single file while iterating
pytest tests/test_perceptivity.py -v

# Coverage
pytest tests/ --cov=perceptome --cov-report=term-missing
```

The test suite includes biological regression checks (paper4.4 / 4.5 / 4.7 a-priori predictions, attractor direction signs, HPA reference v1.4 spot-checks). If a contribution changes one of these numbers, the PR must justify why and reference the empirical work that supports the new value.

## Getting in touch

- Bug reports + feature requests: [GitHub Issues](https://github.com/mool32/perceptome/issues)
- Scientific questions about the framework: open a Discussion or contact tspiro@vaika.org
- Security issues: tspiro@vaika.org (do not file public issue)

By contributing, you agree your contributions are licensed under the same MIT license as the rest of the project.
