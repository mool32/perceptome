"""Build perceptome v0.2.0 module catalog (44 modules: 43 from v0.1 + NPAS4).

Source:
  v0.1 catalog: tool/perceptome/perceptome/data/module_catalog.json (43 modules, internally tagged "version: 0.2")
  NPAS4: Paper 4 Table S4 (neuron-perception 44th candidate, all-criteria pass)

Output:
  tool/perceptome2/perceptome/catalog/data/modules_v03.json (44 modules, version: 0.3)

Adds per-module flags:
  - pan_cellular: bool — does this module operate across all cell types? (44/44 = True for now;
    tissue-specific extensions HNF4A/MYOD/etc. will set False in future v0.4)
  - tissue_bias: list[str] — empty for pan-cellular; for NPAS4 = ["neuron"] as a hint
                  (does NOT exclude module from scoring; scoring should give A_baseline ≈ 0
                  in non-neurons, which is the correct behavior, not an error)
"""

import json
from pathlib import Path

V01_CATALOG = Path("/Users/teo/Desktop/research/perceptual_modules/tool/perceptome/perceptome/data/module_catalog.json")
OUT = Path("/Users/teo/Desktop/research/perceptual_modules/tool/perceptome2/perceptome/catalog/data/modules_v03.json")
OUT.parent.mkdir(parents=True, exist_ok=True)

with open(V01_CATALOG) as f:
    cat = json.load(f)

# Annotate every existing module with pan_cellular flag
for name, mod in cat["modules"].items():
    mod["pan_cellular"] = True
    mod["tissue_bias"] = []

# Add NPAS4 (44th canonical, neuron-perception)
# Source: Paper 4 Table S4. Per Lin et al. 2008 Nature, Spiegel et al. 2014 Cell
# Sensor: Ca2+/CaMKII (synaptic activity); TF: NPAS4 itself (heterodimer with ARNT2)
# Targets: BDNF + 13 neuronal effectors (GABAergic synapse formation, K+ channels)
# Dissociation: HIGH — rapid protein turnover, mRNA peaks transient (Lin 2008)
cat["modules"]["NPAS4"] = {
    "category": "A_exteroceptive",
    "core_genes": ["NPAS4", "ARNT2", "CAMK2A", "CAMK2B", "CALM1"],
    "primary_tf": ["NPAS4"],
    "mii": None,  # to be computed in v0.3 MII recomputation
    "sensor_genes": ["CAMK2A", "CAMK2B", "CALM1", "CALM2", "CALM3"],
    "tf_genes": ["NPAS4", "ARNT2"],
    "feedback_genes": [],  # primary feedback is post-translational (proteasomal turnover),
                            # not transcriptional — kept empty to avoid false signal
    "activity_genes": [
        "BDNF", "NPTX2", "PLK2", "FRMPD4",
        "GABRA1", "GABRA2", "GABRB2", "GABRG2",
        "GPHN", "GAD1", "GAD2",
        "KCNA1", "KCNA2", "NRGN",
    ],
    "dissociation_risk": "HIGH",
    "dissociation_note": "Rapid post-translational turnover (ubiquitin-proteasome). mRNA peaks transient (Lin et al. 2008 Nature). Activity TF targets (BDNF, NPTX2, GABRA1-2, GABRB2, GABRG2, GPHN, GAD1-2, KCNA1-2, NRGN) measure persistent engagement; readiness (NPAS4 mRNA) measures recent induction. ALWAYS report both R and A for neurons.",
    "pan_cellular": True,
    "tissue_bias": ["neuron"],
}

cat["version"] = "0.3"
cat["n_modules"] = len(cat["modules"])
cat["release_notes"] = (
    "v0.3 catalog (perceptome tool v0.2.0): 44 modules. "
    "Adds NPAS4 (neuron-perception 44th canonical, Paper 4 Table S4). "
    "Adds per-module pan_cellular + tissue_bias flags."
)

with open(OUT, "w") as f:
    json.dump(cat, f, indent=2, ensure_ascii=False)

assert cat["n_modules"] == 44
assert "NPAS4" in cat["modules"]
print(f"Wrote {OUT}")
print(f"Modules: {cat['n_modules']}")
print(f"NPAS4 core: {cat['modules']['NPAS4']['core_genes']}")
print(f"NPAS4 activity: {cat['modules']['NPAS4']['activity_genes']}")
