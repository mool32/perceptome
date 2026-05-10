"""Build examples/recipe_experimental_design.ipynb — focused use-case (5-7 cells, no user data)."""

import sys
from pathlib import Path
import nbformat as nbf

OUT = Path(__file__).resolve().parents[1] / "examples" / "recipe_experimental_design.ipynb"


def md(text):
    return nbf.v4.new_markdown_cell(text)


def code(src):
    return nbf.v4.new_code_cell(src)


nb = nbf.v4.new_notebook()
nb.cells = [
    md("""# Recipe — experimental design via capacity-floor predictor

**One question:** I'm planning a perturbation experiment on cell type X to engage pathway Y. Before I run it, will Y actually ramp up, or is it already saturated in cell type X?

**No user data needed** — the bundled HPA reference (154 cell types × 44 modules) answers this directly via the upward-asymmetric capacity-floor predictor.

**Use this for:**
- Drug treatment design — pick the cell line where your target pathway has ramp room, not the one where it's already maxed
- Differentiation / commitment experiments — predict which architecture modules are saturated at the starting cell type
- Avoiding wasted experiments — a pre-experiment check that costs 30 seconds and saves wet-lab time"""),

    code("""import perceptome as pct
import pandas as pd

# Show the bundled HPA reference dimensions
ref = pct.load_hpa_perceptivity()
print(f"HPA reference: {ref['R'].shape[0]} cell types × {ref['R'].shape[1]} modules")
print(f"Sample cell types: {list(ref['R'].index[:5])} ...")
print(f"Sample modules: {list(ref['R'].columns[:5])} ...")"""),

    md("""## Single-cell-type query

**Scenario:** I'm planning a hepatocyte regeneration experiment and want to engage UPR-ATF6, HSF1, and mTOR. Which of these will ramp up?"""),

    code("""pred = pct.predict_engagement(
    starting_cell_type="hepatocytes",
    operation_modules=["UPR-ATF6", "HSF1", "mTOR", "NRF2", "Autophagy"],
)
print(pred[["A_baseline", "headroom", "capacity_floor", "predicted_direction"]].to_string())"""),

    md("""**Reading the result:**

- **`saturated_blocked_up`** = HPA shows this module is already at saturation in this cell type. **Cannot ramp UP** under stimulus. (Active suppression DOWN by specific signaling — e.g. retinoic acid → UPR — IS still possible.)
- **`capacious`** = pathway has ramp room. Magnitude of ramp depends on how strongly your operation engages it (factor 2 of the two-factor framework, not in the tool).
- **`intermediate`** = uncertain — depends on operation intensity.

For the hepatocyte case above: UPR-ATF6 is at A_baseline ≈ 6.3 → blocked up. mTOR similar. **Conclusion**: don't expect strong UP-ramps on these for partial-hepatectomy regeneration; the pathways are already engaged at homeostatic baseline. Pick HSF1 or NRF2 if you need a clear positive signal."""),

    md("""## Cross-cell-type comparison — pick the right model system

**Scenario:** I want to study HSF1 induction. Which cell type gives me the cleanest readout (most ramp room from low baseline)?"""),

    code("""# Loop over a panel of common research cell types
panel = [
    "cardiomyocytes", "hepatocytes", "alveolar cells type 2",
    "fibroblasts", "skeletal myocytes", "brain excitatory neurons",
    "b-cells", "t-cells", "macrophages",
]
rows = []
for ct in panel:
    if ct not in ref["R"].index:
        # Try fuzzy match
        cands = [c for c in ref["R"].index if ct.lower() in c.lower()]
        if not cands:
            continue
        ct = cands[0]
    pred = pct.predict_engagement(ct, ["HSF1"])
    rows.append({
        "cell_type": ct,
        "A_baseline": pred.loc["HSF1", "A_baseline"],
        "headroom": pred.loc["HSF1", "headroom"],
        "capacity_floor": pred.loc["HSF1", "capacity_floor"],
    })
panel_df = pd.DataFrame(rows).sort_values("headroom", ascending=False)
print("HSF1 capacity panel (highest headroom = best ramp candidate):")
print(panel_df.to_string(index=False))"""),

    md("""**Reading the panel:**

- Cells at the top of the headroom column are the best HSF1-engagement candidates — they have low baseline and lots of room to ramp UP under heat shock or proteostatic stress.
- Cells with `saturated_blocked_up` are at the ceiling — picking them means you'll see no signal even with a real biological effect.
- This was used a-priori in Paper 4.4: cardiomyocyte HSF1 was predicted `capacious` (A_baseline ≈ 2.0) → muscle hypertrophy experiment confirmed strong ramp d=+1.29. ✓

## Reference

- The capacity-floor predictor was **closed across paper4.5 + 4.6 + 4.7 + 4.8** (intestinal homeostasis, organoid differentiation, organoid regeneration, organoid retinoid perturbation) as upward-asymmetric.
- Threshold values (`A_baseline > 4.5` = saturated, `< 2.5` = capacious) come from validation across 4 datasets in 3 paradigms.
- For the full scope + paper4.8 PARTIAL refinement (specific signaling can DOWN-suppress saturated modules), see [`docs/SCOPE.md`](../docs/SCOPE.md)."""),
]

nbf.write(nb, OUT)
print(f"wrote {OUT}")
print(f"  cells: {len(nb.cells)}")
