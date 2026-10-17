"""Within-patient deepening pooled across datasets (paper Note S-D).

Pooling the patients of several datasets in one sign test treats them as exchangeable and risks a
Simpson's-paradox artifact. The paper's control is a linear mixed-effects model of the per-patient
malignant-minus-normal difference in core fraction with a dataset random intercept: in its four
low-baseline datasets (two colorectal, kidney, liver; 86 patients) a fixed deepening of +13.0
percentage points (SE 1.87, z = 6.95, one-sided P = 1.9e-12), all four datasets deepening in the mean.
The paper notes the coarseness of a four-group random-effect estimate.

The model is an intercept plus a random intercept per dataset, so it is fitted here directly (REML,
profile over the variance ratio) without a statistics package; ``tests/test_pooled.py`` checks it
against ``statsmodels`` when that is installed.
"""

import numpy as np
import pandas as pd
from scipy import optimize, stats


def _reml_intercept(y, g):
    """REML fit of y_ij = mu + u_i + e_ij, u_i ~ N(0, s_u2), e_ij ~ N(0, s2).

    Returns dict(mu, se, s2, s_u2, lam). With V_i = s2 (I + lam J), the profile REML criterion in
    lam = s_u2 / s2 is (N-1) log Q(lam) + sum_i log(1 + n_i lam) + log sum_i w_i with
    w_i = n_i / (1 + n_i lam) and Q the within sum of squares plus sum_i w_i (ybar_i - mu)^2.
    """
    codes, _ = pd.factorize(g)
    k = int(codes.max()) + 1
    n = np.bincount(codes, minlength=k).astype(float)
    ybar = np.bincount(codes, weights=y, minlength=k) / n
    ssw = float(((y - ybar[codes]) ** 2).sum())
    N = float(n.sum())
    if N < 2:
        raise ValueError("need at least two patients.")

    def parts(lam):
        w = n / (1.0 + n * lam)
        sw = w.sum()
        mu = float((w * ybar).sum() / sw)
        return w, sw, mu, ssw + float((w * (ybar - mu) ** 2).sum())

    def crit(lam):
        _, sw, _, Q = parts(lam)
        return (N - 1) * np.log(max(Q, 1e-300)) + np.log1p(n * lam).sum() + np.log(sw)

    lam = 0.0
    if k > 1:
        grid = np.concatenate([[0.0], np.logspace(-8, 8, 161)])
        vals = np.array([crit(v) for v in grid])
        j = int(vals.argmin())
        lam = float(grid[j])
        if 0 < j < len(grid) - 1:
            lo, hi = np.log(grid[max(j - 1, 1)]), np.log(grid[j + 1])
            res = optimize.minimize_scalar(lambda t: crit(np.exp(t)), bounds=(lo, hi), method="bounded",
                                           options={"xatol": 1e-10})
            if crit(np.exp(res.x)) < crit(lam):
                lam = float(np.exp(res.x))
    _, sw, mu, Q = parts(lam)
    s2 = Q / (N - 1)
    return {"mu": mu, "se": float(np.sqrt(s2 / sw)), "s2": float(s2), "s_u2": float(lam * s2), "lam": lam}


def _long_frame(per_patient, column, dataset):
    if isinstance(per_patient, dict):
        parts = []
        for name, obj in per_patient.items():
            if isinstance(obj, pd.DataFrame) and column not in obj:
                raise ValueError(f"the frame of dataset {name!r} has no column {column!r}.")
            s = obj[column] if isinstance(obj, pd.DataFrame) else pd.Series(obj)
            parts.append(pd.DataFrame({"dataset": name, "delta": s.astype(float).values}))
        long = pd.concat(parts, ignore_index=True)
    else:
        df = pd.DataFrame(per_patient)
        if dataset is None or dataset not in df or column not in df:
            raise ValueError(f"a DataFrame needs the columns {column!r} and the dataset column (dataset=...).")
        long = pd.DataFrame({"dataset": df[dataset].values, "delta": df[column].astype(float).values})
    return long.dropna()


def pooled_within_patient(per_patient, column="delta", dataset=None):
    """Mixed-effects estimate of the within-patient deepening pooled over datasets (Note S-D).

    Fits the per-patient malignant-minus-normal difference as an intercept plus a dataset random
    intercept (REML) and tests the intercept with a Wald z, as the paper does.

    Parameters
    ----------
    per_patient : dict | DataFrame
        Either ``{dataset name: per_patient frame of within_patient_occupancy (or a Series of
        differences)}``, or one DataFrame with a row per patient, the difference in ``column`` and
        the dataset label in the column named by ``dataset``.
    column : str
        The column holding the difference (``within_patient_occupancy(...)['per_patient']['delta']``,
        a fraction; use ``delta_matched`` for the complexity-matched difference).
    dataset : str | None
        Name of the dataset column when ``per_patient`` is a DataFrame.

    Returns
    -------
    dict
        estimate, se             fixed deepening and its standard error, as fractions
        estimate_points, se_points   the same in percentage points (the paper's unit)
        z, p_one_sided, p_two_sided  Wald z = estimate / se; P of a deepening (one-sided) as in the paper
        sigma2_dataset, sigma2_residual   estimated variance components (dataset, patient)
        n_patients, n_datasets
        per_dataset              DataFrame: patients, patients deeper, mean difference (points)

    Notes
    -----
    With few datasets the between-dataset variance is estimated coarsely (the paper's four are two
    colorectal cohorts, kidney and liver), so the result is carried by the datasets that have
    patients to spare; the one-sided P does not propagate the uncertainty of the variance estimate.
    """
    long = _long_frame(per_patient, column, dataset)
    if len(long) < 2:
        raise ValueError("need at least two patients with a difference.")
    fit = _reml_intercept(long["delta"].values.astype(float), long["dataset"].values)
    z = fit["mu"] / fit["se"] if fit["se"] > 0 else float(np.sign(fit["mu"]) * np.inf)
    per = long.groupby("dataset", sort=True)["delta"].agg(
        n_patients="size", n_deeper=lambda d: int((d > 0).sum()), mean_delta=lambda d: 100 * float(d.mean()))
    per = per.rename(columns={"mean_delta": "mean_delta_points"})
    return {
        "estimate": fit["mu"], "se": fit["se"],
        "estimate_points": 100 * fit["mu"], "se_points": 100 * fit["se"],
        "z": float(z), "p_one_sided": float(stats.norm.sf(z)), "p_two_sided": float(2 * stats.norm.sf(abs(z))),
        "sigma2_dataset": fit["s_u2"], "sigma2_residual": fit["s2"],
        "n_patients": int(len(long)), "n_datasets": int(long["dataset"].nunique()),
        "per_dataset": per,
    }
