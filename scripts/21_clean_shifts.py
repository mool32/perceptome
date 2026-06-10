#!/usr/bin/env python3
"""Phase 2b/#1 — calibration-independent tumor−normal SHIFT vectors.

A within-dataset paired shift (malignant centroid − matched-normal centroid, same
cell type, same scoring) cancels the per-cell scale offset, so its projection onto
the main-sequence axis is calibration-robust (unlike absolute position). Computed
for every cancer with a within-dataset normal counterpart available locally:
HCC (Sun2021: Hepatocyte Tumor vs Normal) and BRCA (Wu2021: Cancer vs Normal
Epithelial, overall + per subtype).
"""
import json
import sys
import warnings
from pathlib import Path

import numpy as np
import scanpy as sc

warnings.filterwarnings("ignore")
ROOT = Path("/Users/teo/Desktop/research/perceptual_modules")
TOOL = ROOT / "tool" / "perceptome2"
sys.path.insert(0, str(TOOL))
import perceptome as pct  # noqa: E402
from perceptome.eigenspace.project import _load_reference  # noqa: E402

ref = _load_reference(None); mo = ref["module_order"]; npc = ref["n_pcs"]
std_ref = np.array(ref["std_per_module"]); std_ref[std_ref < 1e-10] = 1.0
loadings = np.array([[ref["loadings"][f"PC{j+1}"].get(m, 0.0) for j in range(npc)] for m in mo])

cc = json.load(open(ROOT / "transition_vectors/results/per_pc_decomposition.json"))["class_centroids"]
Mtr = np.array([cc[n] for n in cc]); Mtr = Mtr / np.linalg.norm(Mtr, axis=1, keepdims=True)
ms = np.linalg.svd(Mtr, full_matrices=False)[2][0]
if (Mtr[list(cc).index("TUMOR_TME")] @ ms) < 0:
    ms = -ms

def cos(a, b): return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-12))

_core = set(pct.get_genes("HSF1", "core") + pct.get_genes("p53", "core"))

def score(adata):
    if "feature_name" in adata.var.columns and len(_core & set(adata.var_names)) == 0:
        adata.var_names = adata.var["feature_name"].astype(str).values
        adata.var_names_make_unique()
    if float(adata.X[:200].max()) > 50:
        sc.pp.normalize_total(adata, target_sum=1e4); sc.pp.log1p(adata)
    return pct.score_modules(adata)["scores"].reindex(columns=mo)

def shift_cos(tum_scores, nor_scores):
    d = tum_scores.mean().values - nor_scores.mean().values     # 44-module shift
    v9 = (d / std_ref) @ loadings
    return cos(v9, ms), float(np.linalg.norm(d))

results = []

# ---- HCC (Sun2021) ----
a = sc.read_h5ad(ROOT / "results/cancer_sun2021_hcc/sun2021.h5ad")
s = score(a)
hep = a.obs["celltype"].values == "Hepatocyte"
tum = s[hep & (a.obs["site"].values == "Tumor")]
nor = s[hep & (a.obs["site"].values == "Normal")]
c, mag = shift_cos(tum, nor)
results.append(("HCC (Hepatocyte)", len(tum), len(nor), c, mag))
del a, s

# ---- BRCA (Wu2021): Cancer vs Normal Epithelial, overall + per subtype ----
a = sc.read_h5ad(ROOT / "results/cancer_validation/wu2021_brca.h5ad")
s = score(a)
maj = a.obs["celltype_major"].astype(str).values
sub = a.obs["subtype"].astype(str).values
nor_all = s[maj == "Normal Epithelial"]
can_all = s[maj == "Cancer Epithelial"]
c, mag = shift_cos(can_all, nor_all)
results.append(("BRCA (Cancer Epi, all)", len(can_all), len(nor_all), c, mag))
for st in ["ER+", "HER2+", "TNBC"]:
    can = s[(maj == "Cancer Epithelial") & (sub == st)]
    if len(can) >= 30:
        c, mag = shift_cos(can, nor_all)
        results.append((f"BRCA {st}", len(can), len(nor_all), c, mag))
del a, s

# ---- BCC (Census bcc_full): basal cells, carcinoma vs normal ----
a = sc.read_h5ad(ROOT / "paper6_immune_3layer/data/bcc_full.h5ad")
s = score(a)
basal = a.obs["cell_type"].astype(str).isin(["basal cell", "basal cell of epidermis"]).values
dis = a.obs["disease"].astype(str).values
tum = s[basal & (dis == "basal cell carcinoma")]
nor = s[basal & (dis == "normal")]
c, mag = shift_cos(tum, nor)
results.append(("BCC (basal cell)", len(tum), len(nor), c, mag))
del a, s

# ---- report ----
print("Calibration-independent tumor−normal shift on the main-sequence axis")
print("(cos>0 = malignant shift toward the chronic-reorganization pole)\n")
print(f"  {'cancer':26s} {'n_tum':>7s} {'n_norm':>7s}  cos(shift, main-seq)   ‖shift‖")
for name, nt, nn, c, mag in results:
    print(f"  {name:26s} {nt:>7d} {nn:>7d}        {c:+.2f}            {mag:.3f}")

print("\nreference: aggregate TUMOR_TME transition = +0.85; all clean shifts should be > 0 "
      "(chronic-ward),\n           independent of per-cell calibration (the offset cancels in tumor−normal).")
