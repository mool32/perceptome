#!/usr/bin/env python3
"""Phase 2b — geometric characterization of the perceptome manifold.

Three questions toward a unified, interpretable coordinate system (the cellular
Hertzsprung–Russell analog): is the 154-cell manifold (1) LOW-DIMENSIONAL,
(2) CURVED, and (3) do distinct biological processes live as trajectories on the
SAME manifold (does cancer convergence ≈ reverse differentiation)?

Uses the single_cell_scaled calibration (the native space of the cancer-
convergence finding, paper3 v7) and the bundled 9-PC eigenspace loadings.
"""
import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
import perceptome as pct  # noqa: E402
from perceptome.eigenspace.project import _load_reference  # noqa: E402
from perceptome.eigenspace.manifold import DiffusionMap, _sqdist  # noqa: E402

# ---------- load manifold (sc-scaled, projected through 9-PC loadings) ----------
ref = _load_reference(None)
mo = ref["module_order"]
npc = ref["n_pcs"]
mean = np.array(ref["mean_per_module"]); std = np.array(ref["std_per_module"])
std = np.where(std < 1e-10, 1.0, std)
loadings = np.array([[ref["loadings"][f"PC{j+1}"].get(m, 0.0) for j in range(npc)] for m in mo])
eigvals = np.array(ref["eigenvalues"])

# Canonical interpretable space = pseudobulk readiness z-scored against the SAME
# per-module mean/std the eigenspace was built from (reference_v03). (single_cell_scaled
# is a per-cell-projection calibration for tumor data; z-scoring it against pseudobulk
# mean/std collapses everything to the erythrocyte pole — the v0.2.3 artifact.)
hpa = pct.load_hpa_perceptivity(mode="pseudobulk")
R = hpa["R"].reindex(columns=mo)
cells = list(R.index)
Xz = (R.values.astype(float) - mean) / std       # z-scored module space (154 × 44), self-consistent
PCA = Xz @ loadings                                # 154 × 9 canonical coordinates
N = len(cells)

print("=" * 78)
print("1. DIMENSIONALITY  — is the manifold low-dimensional?")
print("=" * 78)
# (a) linear: participation ratio of PCA eigenvalues
pr = (eigvals.sum() ** 2) / (eigvals ** 2).sum()
ev9 = eigvals[:npc]
print(f"  Kaiser linear dims (eigenvalue>1):      {npc}")
print(f"  participation ratio (all 44 eigvals):   {pr:.2f}  effective linear dims")
print(f"  variance in PC1 / PC1-3 / PC1-9:         {ev9[0]/eigvals.sum():.1%} / "
      f"{ev9[:3].sum()/eigvals.sum():.1%} / {ev9.sum()/eigvals.sum():.1%}")

# (b) intrinsic dim, TwoNN (Facco 2017) on the z-scored module space
def twonn(X):
    from scipy.spatial import cKDTree
    d, _ = cKDTree(X).query(X, k=3)
    mu = d[:, 2] / np.maximum(d[:, 1], 1e-12)
    mu = np.sort(mu[np.isfinite(mu) & (mu > 1)])
    n = len(mu)
    F = np.arange(1, n + 1) / n
    keep = F < 0.9
    x = np.log(mu[keep]); y = -np.log(1 - F[keep] + 1e-12)
    return float(np.sum(x * y) / np.sum(x * x))

print(f"  TwoNN intrinsic dim (module space):      {twonn(Xz):.2f}")
print(f"  TwoNN intrinsic dim (9-PC space):        {twonn(PCA):.2f}")

# (c) diffusion eigenvalue spectrum (nonlinear) — eigengap
dm = DiffusionMap(n_components=20).fit(Xz)
dvals = dm.evals_[1:21]
gaps = dvals[:-1] - dvals[1:]
print(f"  diffusion eigenvalues λ1..λ6:            {np.round(dvals[:6], 3)}")
print(f"  largest eigengap after component:        {int(np.argmax(gaps)) + 1}  "
      f"(gap={gaps.max():.3f})")

print("\n" + "=" * 78)
print("2. CURVATURE  — is the manifold flat or curved?")
print("=" * 78)
from sklearn.neighbors import kneighbors_graph
from scipy.sparse.csgraph import shortest_path
from scipy.stats import pearsonr

def curvature(X, k):
    g = kneighbors_graph(X, n_neighbors=k, mode="distance")
    geo = shortest_path(g, method="D", directed=False)
    if not np.isfinite(geo).all():
        return None
    iu = np.triu_indices(len(X), 1)
    eucl = np.sqrt(_sqdist(X))[iu]; gg = geo[iu]
    far = eucl > np.percentile(eucl, 75)
    return pearsonr(eucl, gg)[0] ** 2, (gg[far] / eucl[far]).mean()

# null control: flat-random clouds matched in N and effective dimension, to calibrate
# how much of the geo/eucl ratio is just sparse-sampling vs genuine curvature.
rng = np.random.default_rng(0)
flat_iso = rng.standard_normal((N, int(round(pr))))                    # isotropic flat
shuffled = np.column_stack([rng.permutation(Xz[:, j]) for j in range(Xz.shape[1])])  # destroys manifold, keeps margins
for k in (8, 12):
    real = curvature(Xz, k)
    fi = curvature(flat_iso, k)
    sh = curvature(shuffled, k)
    if real is None:
        print(f"  k={k}: kNN graph disconnected — skip"); continue
    print(f"  k={k}:  REAL   R²={real[0]:.3f}  ratio={real[1]:.2f}   |  "
          f"flat-random ratio={fi[1]:.2f}   shuffled-features ratio={sh[1]:.2f}")
print("  curvature is REAL only if real ratio > flat/shuffled baselines "
      "(else it's sparse-sampling).")

print("\n" + "=" * 78)
print("3. PRINCIPAL AXIS  — is there a 'main sequence' backbone?")
print("=" * 78)
pc1 = PCA[:, 0]
order = np.argsort(pc1)
print(f"  PC1 carries {ev9[0]/eigvals.sum():.1%} of variance (dominant axis)")
print("  PC1 low pole (5):  ", [cells[i] for i in order[:5]])
print("  PC1 high pole (5): ", [cells[i] for i in order[-5:]])

print("\n" + "=" * 78)
print("4. TRAJECTORIES ON ONE MANIFOLD  — cancer convergence vs differentiation")
print("=" * 78)
def centroid(names):
    idx = [cells.index(c) for c in names if c in cells]
    return PCA[idx].mean(0), idx

# v7 convergence core + neighborhood
duet = ["megakaryocytes", "cytotrophoblasts"]
embryonic = ["cytotrophoblasts", "migrating cytotrophoblasts", "syncytiotrophoblasts",
             "extravillous trophoblasts", "oocytes", "early primary spermatocytes",
             "late primary spermatocytes", "early spermatids", "late spermatids",
             "differentiating spermatogonia", "undifferentiated spermatogonia"]
# terminal differentiated pole (low perception-breadth, post-mitotic specialists)
terminal = ["erythrocytes", "cardiomyocytes", "hepatocytes", "brain excitatory neurons",
            "enterocytes", "adipocytes", "rod photoreceptor cells"]

c_duet, _ = centroid(duet)
c_emb, _ = centroid(embryonic)
c_term, _ = centroid(terminal)

# differentiation axis: embryonic -> terminal ; cancer convergence: toward duet/embryonic
diff_axis = c_term - c_emb
cancer_dir = np.array(pct.load_attractor_direction()["attractor_direction_eigenspace"])[:npc]

def cos(a, b):
    return float(a @ b / (np.linalg.norm(a) * np.linalg.norm(b) + 1e-12))

print(f"  duet mutual cosine (9-PC eigenspace):    "
      f"{cos(PCA[cells.index('megakaryocytes')], PCA[cells.index('cytotrophoblasts')]):.2f}  "
      f"(paper v7 reports ~0.70 for the convergence neighborhood)")
print(f"  embryonic-pole vs terminal-pole on PC1:  {c_emb[0]:.2f}  vs  {c_term[0]:.2f}")
print(f"  cos(cancer_direction, differentiation_axis):     {cos(cancer_dir, diff_axis):+.2f}")
print(f"  cos(cancer_direction, REVERSE differentiation):  {cos(cancer_dir, -diff_axis):+.2f}")
print("  → if cancer aligns with REVERSE differentiation, de-diff retraces the path back")

# how low-dimensional are the process directions? share of each direction in PC1..PC4
for nm, v in [("differentiation_axis", diff_axis), ("cancer_direction", cancer_dir)]:
    w = v / (np.linalg.norm(v) + 1e-12)
    print(f"  {nm:22s} energy in PC1..PC4: {np.round(w[:4]**2, 2)}  "
          f"(sum={ (w[:4]**2).sum():.2f} of 1.0)")

# embryonic co-localization with duet (reproduce paper's 9/11 in top-30 as geometry check)
dcent, _ = centroid(duet)
cos_to_duet = np.array([cos(PCA[i], dcent) for i in range(N)])
rank = np.argsort(-cos_to_duet)
emb_ranks = sorted(int(np.where(rank == cells.index(c))[0][0]) + 1 for c in embryonic if c in cells)
print(f"\n  embryonic cells' rank (of 154) by cosine to duet centroid: {emb_ranks}")
print(f"  embryonic in top-30: {sum(r<=30 for r in emb_ranks)}/{len(emb_ranks)} "
      f"(paper v7: 9/11, Fisher P=5.6e-6)")
