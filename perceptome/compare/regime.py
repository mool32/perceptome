"""4-quadrant infrastructure regime classification.

Rule (paper3 §3.5):
  perception_delta > 0, infrastructure_delta > 0  →  supply_chain    (healthy scaling)
  perception_delta < 0, infrastructure_delta > 0  →  firefighting     (damage response)
  perception_delta < 0, infrastructure_delta < 0  →  collapse         (system failure)
  perception_delta > 0, infrastructure_delta < 0  →  unsupported      (unstable activation)

Module assignments (extended from v0.1.0 to include NPAS4 in perception):
  INFRASTRUCTURE: HSF1, UPR-IRE1, UPR-PERK, UPR-ATF6, NRF2, mTOR, Autophagy
  PERCEPTION:    NF-κB, ERK/MAPK, JAK-STAT, cAMP/CREB, NFAT, Wnt, Notch, Hippo,
                 TGF-β, p53, Calcium, cGAS-STING, Type I IFN, BMP, Hedgehog,
                 PI3K/PTEN, NPAS4
"""

import pandas as pd

_INFRA_MODULES = (
    "HSF1", "UPR-IRE1", "UPR-PERK", "UPR-ATF6", "NRF2", "mTOR", "Autophagy",
)
_PERCEPTION_MODULES = (
    "NF-κB", "ERK/MAPK", "JAK-STAT", "cAMP/CREB", "NFAT",
    "Wnt", "Notch", "Hippo", "TGF-β", "p53", "Calcium",
    "cGAS-STING", "Type I IFN", "BMP", "Hedgehog", "PI3K/PTEN",
    "NPAS4",
)


def infrastructure_regime(delta_modules):
    """Classify infrastructure regime from per-module deltas.

    Parameters
    ----------
    delta_modules : DataFrame | Series
        Module deltas (e.g., disease − control). If DataFrame, uses 'delta' column.

    Returns
    -------
    dict
        regime               supply_chain | firefighting | collapse | unsupported | neutral
        perception_delta     mean delta across PERCEPTION modules
        infrastructure_delta mean delta across INFRASTRUCTURE modules
        description          one-line explanation
    """
    if isinstance(delta_modules, pd.DataFrame):
        deltas = delta_modules["delta"] if "delta" in delta_modules.columns else delta_modules.iloc[:, 0]
    else:
        deltas = delta_modules

    perc_avail = [m for m in _PERCEPTION_MODULES if m in deltas.index]
    infra_avail = [m for m in _INFRA_MODULES if m in deltas.index]

    perc_delta = float(deltas[perc_avail].mean()) if perc_avail else 0.0
    infra_delta = float(deltas[infra_avail].mean()) if infra_avail else 0.0

    if perc_delta > 0 and infra_delta > 0:
        regime, desc = "supply_chain", "Healthy scaling: perception and infrastructure both elevated"
    elif perc_delta < 0 and infra_delta > 0:
        regime, desc = "firefighting", "Damage response: infrastructure fighting fires while perception declines"
    elif perc_delta < 0 and infra_delta < 0:
        regime, desc = "collapse", "System failure: both perception and infrastructure declining"
    elif perc_delta > 0 and infra_delta < 0:
        regime, desc = "unsupported", "Unsupported activation: perception up without infrastructure support (unstable)"
    else:
        regime, desc = "neutral", "Minimal change"

    return {
        "regime": regime,
        "perception_delta": round(perc_delta, 6),
        "infrastructure_delta": round(infra_delta, 6),
        "description": desc,
    }
