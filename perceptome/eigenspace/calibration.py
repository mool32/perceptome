"""Reference calibrations — how module scores are standardized before projection.

The reference coordinate system is built once from the 154 normal HPA cell types
and never refit. Placing a cell in it takes three steps (Spiro 2026, Methods and
Note S1): score the 44 modules (``mean_raw``), standardize every module against
the reference, project through the fixed 44 x 9 loading matrix.

Step two depends on the *scale* of the reference, and the paper documents two
calibrations of the same gene panel:

  pseudobulk           reference scored on population-averaged HPA expression
                       (mean log1p(nCPM)). The eigenspace is built on it and it is
                       the calibration bundled with this package.
  single_cell_scaled   reference re-derived under the same
                       normalize_total(1e4) + log1p preprocessing as the per-cell
                       data, so per-cell z-scores are commensurate. Required for
                       per-cell projection: under the pseudobulk calibration tumour
                       cells collapse onto low-complexity targets (Note S-B), and so
                       do independent normal single cells (99.6% onto erythroid
                       targets, Note S-H).

A :class:`Calibration` bundles the module means / SDs, the loading matrix and
(optionally) the reference score matrix, from which the 154 normal-cell centroids
and the shape-only variant (Note S-K iv) are derived.
"""

import json
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

_DATA_DIR = Path(__file__).parent / "data"
_SINGLE_CELL_FILE = _DATA_DIR / "calibration_single_cell_scaled_v03.json"

PSEUDOBULK = "pseudobulk"
SINGLE_CELL_SCALED = "single_cell_scaled"

# Median standardized module score (in reference SDs) below which input looks like
# per-cell scores projected with a pseudobulk calibration.
_SCALE_WARN_Z = -1.5


class Calibration:
    """Standardization statistics + loadings (+ reference scores) for one scaling.

    Parameters
    ----------
    name : str
    module_order : sequence of str
    mean, std : array-like (n_modules,)
        Reference mean / SD per module (population SD, ddof=0, as in ``rebuild``).
    loadings : DataFrame (n_modules x n_pcs)
        Fixed eigenvector matrix; index = module names, columns = PC1..PCk.
    reference_scores : DataFrame | None
        154 x 44 reference matrix under this scaling. Needed for ``centroids`` and
        ``shape_only_calibration``.
    shape_only : bool
        If True, every input row (and the reference) has its own mean across the
        44 modules subtracted before standardization (Note S-K iv).
    """

    def __init__(self, name, module_order, mean, std, loadings,
                 reference_scores=None, shape_only=False):
        self.name = name
        self.module_order = list(module_order)
        self.mean = pd.Series(np.asarray(mean, dtype=float), index=self.module_order)
        std = np.asarray(std, dtype=float).copy()
        std[std < 1e-10] = 1.0
        self.std = pd.Series(std, index=self.module_order)
        self.loadings = loadings.loc[self.module_order].astype(float)
        self.reference_scores = reference_scores
        self.shape_only = bool(shape_only)
        self._centroids = None

    # -- construction -------------------------------------------------------

    @classmethod
    def from_reference(cls, ref, name=PSEUDOBULK, reference_scores=None):
        """Build from an eigenspace reference dict (``reference_v03.json`` layout)."""
        if "mean_per_module" not in ref or "std_per_module" not in ref:
            raise ValueError(
                "Eigenspace reference carries no calibration statistics "
                "(mean_per_module / std_per_module); rebuild it with pct.rebuild() "
                "or project with standardize=False."
            )
        order = list(ref["module_order"])
        return cls(name, order, ref["mean_per_module"], ref["std_per_module"],
                   _loadings_frame(ref, order), reference_scores=reference_scores)

    @classmethod
    def from_scores(cls, reference_scores, loadings, name):
        """Build from a cell-type x module reference matrix (statistics derived from it)."""
        order = list(loadings.index)
        ref = reference_scores[order].astype(float)
        return cls(name, order, ref.mean(axis=0).values, ref.std(axis=0, ddof=0).values,
                   loadings, reference_scores=ref)

    # -- core operations ----------------------------------------------------

    @property
    def n_pcs(self):
        return self.loadings.shape[1]

    def standardize(self, scores):
        """Module scores -> reference SD units.

        Modules absent from ``scores`` are placed at the reference mean (z = 0),
        i.e. they neither pull nor push a projection.
        """
        present = [m for m in self.module_order if m in scores.columns]
        X = scores[present].astype(float)
        if self.shape_only:
            X = X.sub(X.mean(axis=1), axis=0)
        Z = (X - self.mean[present]) / self.std[present]
        return Z.reindex(columns=self.module_order, fill_value=0.0)

    def project(self, scores, n_pcs=None, missing_modules=None):
        """Standardize ``scores`` and project through the loadings -> samples x PCs."""
        k = self.n_pcs if n_pcs is None else min(int(n_pcs), self.n_pcs)
        Z = self.standardize(scores)
        if missing_modules:
            Z.loc[:, [m for m in missing_modules if m in Z.columns]] = 0.0
        coords = Z.values @ self.loadings.values[:, :k]
        return pd.DataFrame(coords, index=scores.index,
                            columns=list(self.loadings.columns[:k]))

    def project_shift(self, delta, n_pcs=None):
        """Project module-score *differences* (scale by 1/SD, no centering).

        ``project(a) - project(b) == project_shift(a - b)`` for rows ``a``, ``b``.
        """
        k = self.n_pcs if n_pcs is None else min(int(n_pcs), self.n_pcs)
        single = isinstance(delta, pd.Series)
        D = delta.to_frame().T if single else delta
        present = [m for m in self.module_order if m in D.columns]
        Z = (D[present].astype(float) / self.std[present]).reindex(
            columns=self.module_order, fill_value=0.0)
        out = pd.DataFrame(Z.values @ self.loadings.values[:, :k], index=D.index,
                           columns=list(self.loadings.columns[:k]))
        return out.iloc[0].rename("shift") if single else out

    @property
    def centroids(self):
        """Normal-cell centroids in the eigenspace (reference cell types x PCs)."""
        if self.reference_scores is None:
            raise ValueError(
                f"Calibration {self.name!r} was built from statistics only and holds no "
                "reference score matrix, so it has no cell-type centroids."
            )
        if self._centroids is None:
            self._centroids = self.project(self.reference_scores)
        return self._centroids

    def shape_only_calibration(self):
        """Shape-only variant: subtract each cell's own mean across the 44 modules.

        The reference is centred the same way and its mean / SD recomputed on the
        centred matrix; loadings, projection and assignment are otherwise unchanged
        (Note S-K iv).
        """
        if self.reference_scores is None:
            raise ValueError("shape_only needs the reference score matrix.")
        if self.shape_only:
            return self
        ref = self.reference_scores[self.module_order].astype(float)
        centred = ref.sub(ref.mean(axis=1), axis=0)
        return Calibration(
            f"{self.name}+shape_only", self.module_order,
            centred.mean(axis=0).values, centred.std(axis=0, ddof=0).values,
            self.loadings, reference_scores=centred, shape_only=True,
        )

    # -- persistence --------------------------------------------------------

    def to_dict(self):
        ref = self.reference_scores
        return {
            "name": self.name,
            "shape_only": self.shape_only,
            "module_order": self.module_order,
            "mean": self.mean.tolist(),
            "std": self.std.tolist(),
            "loadings": {pc: self.loadings[pc].to_dict() for pc in self.loadings.columns},
            "reference_scores": None if ref is None else {
                "index": list(ref.index), "columns": list(ref.columns),
                "values": ref.values.astype(float).tolist(),
            },
        }

    @classmethod
    def from_dict(cls, d):
        order = list(d["module_order"])
        loadings = pd.DataFrame({pc: pd.Series(v) for pc, v in d["loadings"].items()}).loc[order]
        rs = d.get("reference_scores")
        ref = None if rs is None else pd.DataFrame(rs["values"], index=rs["index"], columns=rs["columns"])
        return cls(d["name"], order, d["mean"], d["std"], loadings,
                   reference_scores=ref, shape_only=d.get("shape_only", False))

    def save(self, path):
        with open(path, "w") as f:
            json.dump(self.to_dict(), f, ensure_ascii=False)

    def __repr__(self):
        ref = "" if self.reference_scores is None else f", {len(self.reference_scores)} reference cell types"
        return f"Calibration({self.name!r}, {len(self.module_order)} modules, {self.n_pcs} PCs{ref})"


def _loadings_frame(ref, module_order):
    n = ref["n_pcs"]
    cols = {f"PC{j + 1}": [ref["loadings"][f"PC{j + 1}"].get(m, 0.0) for m in module_order]
            for j in range(n)}
    return pd.DataFrame(cols, index=module_order)


def load_calibration(name=PSEUDOBULK, reference=None):
    """Load a reference calibration.

    Parameters
    ----------
    name : 'pseudobulk' | 'single_cell_scaled' | path to a calibration JSON
        ``pseudobulk`` is the bundled calibration (eigenspace reference + the
        154 x 44 HPA readiness matrix). ``single_cell_scaled`` is the calibration
        the paper uses for per-cell data; it is loaded from
        ``eigenspace/data/calibration_single_cell_scaled_v03.json`` and has to be
        built from the HPA single-cell-type table first
        (``scripts/13_build_single_cell_calibration.py``).
    reference : dict | str | None
        Eigenspace reference (loadings) to pair with ``pseudobulk``.
    """
    from .project import _load_reference

    if isinstance(name, Calibration):
        return name
    if name == PSEUDOBULK:
        from ..perceptivity.reference import load_hpa_perceptivity
        ref = _load_reference(reference) if (reference is None or isinstance(reference, str)) else reference
        order = list(ref["module_order"])
        scores = load_hpa_perceptivity()["R"][order].astype(float)
        cal = Calibration.from_reference(ref, PSEUDOBULK, reference_scores=scores)
        return cal
    if name == SINGLE_CELL_SCALED:
        if not _SINGLE_CELL_FILE.exists():
            raise FileNotFoundError(
                "The single_cell_scaled calibration is not bundled with this package: it has to be "
                "derived from the Human Protein Atlas single-cell-type table "
                "(rna_single_cell_type.tsv), which is not redistributed here. Build it with "
                "scripts/13_build_single_cell_calibration.py <path-to-HPA-tsv> and check it against "
                "Table S7 of the paper (pct.convergence.check_calibration_against_table_s7). "
                "See docs/PAPER_MAP.md."
            )
        name = _SINGLE_CELL_FILE
    path = Path(name)
    if not path.exists():
        raise FileNotFoundError(f"Calibration file not found: {path}")
    with open(path) as f:
        return Calibration.from_dict(json.load(f))


def hpa_reference_scores(expression, catalog=None, gene_set="core",
                         scaling=PSEUDOBULK, target_sum=1e4, library_size=None):
    """Cell-type x module reference matrix under a chosen scaling.

    Parameters
    ----------
    expression : DataFrame (genes x cell types)
        HPA single-cell-type expression (nCPM or read counts), indexed by gene symbol.
    scaling : 'pseudobulk' | 'single_cell_scaled'
        ``pseudobulk``: mean over a module's genes of ``log1p(expression)`` — what the
        bundled 154 x 44 matrix is (mean log1p(nCPM)).
        ``single_cell_scaled``: each cell type is first scaled to ``target_sum`` total
        (``normalize_total(1e4)``), then ``log1p`` — the preprocessing applied to the
        per-cell data (Note S1).
    library_size : Series | None
        Per-cell-type total over *all* genes. Defaults to ``expression.sum()``, which is
        only correct if ``expression`` holds the full gene set.
    """
    from ..catalog import load_catalog, get_genes

    cat = load_catalog() if catalog is None else catalog
    expr = expression.astype(float)
    if scaling == SINGLE_CELL_SCALED:
        if library_size is None:
            if expr.shape[0] < 5000:
                warnings.warn(
                    f"expression holds only {expr.shape[0]} genes; normalize_total needs the total over "
                    "all genes. Pass library_size= or the full gene table.", UserWarning, stacklevel=2)
            library_size = expr.sum(axis=0)
        expr = expr.mul(target_sum / pd.Series(library_size).reindex(expr.columns), axis=1)
    elif scaling != PSEUDOBULK:
        raise ValueError(f"scaling must be 'pseudobulk' or 'single_cell_scaled', got {scaling!r}")
    logged = np.log1p(expr)

    modules = sorted(cat["modules"])
    out = {}
    for m in modules:
        genes = [g for g in get_genes(m, gene_set, cat) if g in logged.index]
        out[m] = logged.loc[genes].mean(axis=0) if genes else pd.Series(np.nan, index=logged.columns)
    return pd.DataFrame(out)[modules]


def read_hpa_single_cell_type(path, genes=None, value_column="nCPM"):
    """Stream the HPA ``rna_single_cell_type.tsv`` table.

    Returns ``(expression, library_size)``: ``expression`` is genes x cell types for the requested
    ``genes`` (all genes if None); ``library_size`` is the per-cell-type total over *all* genes in
    the file, which ``normalize_total`` needs. Where several rows share a gene symbol the last one
    wins, as in ``scripts/02_build_hpa_perceptivity.py``.

    Columns are found by header name (``Gene name``, ``Cell type``, ``value_column``); a four-column
    file without those names is read positionally (gene, gene name, cell type, nCPM).
    """
    keep = None if genes is None else set(genes)
    expr, totals = {}, {}
    with open(path, encoding="utf-8") as f:
        header = f.readline().rstrip("\n").split("\t")
        col = {name: i for i, name in enumerate(header)}
        if {"Gene name", "Cell type", value_column} <= set(col):
            i_gene, i_ct, i_val = col["Gene name"], col["Cell type"], col[value_column]
        elif len(header) == 4:
            i_gene, i_ct, i_val = 1, 2, 3
        else:
            raise ValueError(f"Expected columns 'Gene name', 'Cell type' and {value_column!r}; found {header}")
        need = max(i_gene, i_ct, i_val)
        for line in f:
            parts = line.rstrip("\n").split("\t")
            if len(parts) <= need:
                continue
            try:
                v = float(parts[i_val])
            except ValueError:
                continue
            ct, g = parts[i_ct], parts[i_gene]
            totals[ct] = totals.get(ct, 0.0) + v
            if keep is None or g in keep:
                expr.setdefault(g, {})[ct] = v
    cell_types = sorted(totals)
    expression = pd.DataFrame(expr, index=cell_types).T
    return expression, pd.Series(totals).reindex(cell_types)


def scale_diagnostic(z):
    """Median over samples of the mean standardized score across modules (reference SDs)."""
    return float(np.nanmedian(np.nanmean(np.asarray(z, dtype=float), axis=1)))


def warn_if_scale_mismatch(z, calibration, stacklevel=3):
    """Warn when scores sit far below a pseudobulk reference (looks like per-cell input)."""
    if getattr(calibration, "name", "") != PSEUDOBULK or len(z) == 0:
        return None
    d = scale_diagnostic(z)
    if d < _SCALE_WARN_Z:
        warnings.warn(
            f"Module scores sit {abs(d):.1f} reference SDs below the pseudobulk calibration on average. "
            "That is the signature of per-cell scores (normalize_total(1e4)+log1p) projected with the "
            "population-averaged reference, which collapses cells onto low-complexity reference types "
            "(paper Notes S-B, S-H). For per-cell data use a single_cell_scaled calibration "
            "(pct.load_calibration('single_cell_scaled')). Ignore this if the scores are "
            "pseudobulk-scale aggregates.", UserWarning, stacklevel=stacklevel)
    return d
