"""PASS/FAIL/ARTIFACT scorecard for perturbation analyses."""

from dataclasses import dataclass, field, asdict
from typing import List

from .nulls import (
    random_200_panel, housekeeping_panel, cell_cycle_panel,
    log2fc_perturbation_vs_control, _load_panels,
)


@dataclass
class CheckResult:
    name: str
    value: float
    threshold: float
    direction: str            # 'lt' (must be < threshold) or 'gt' (must be > threshold)
    verdict: str              # PASS | FAIL | ARTIFACT | NA
    notes: str = ""


@dataclass
class Scorecard:
    overall_verdict: str
    checks: List[CheckResult] = field(default_factory=list)

    def to_dict(self):
        return {
            "overall_verdict": self.overall_verdict,
            "checks": [asdict(c) for c in self.checks],
        }


def _classify(value, threshold, direction):
    if value != value:  # NaN
        return "NA"
    if direction == "lt":
        return "PASS" if abs(value) < threshold else "FAIL"
    if direction == "gt":
        return "PASS" if abs(value) > threshold else "FAIL"
    return "NA"


def validate_perturbation(
    adata,
    condition_col,
    perturbation_value,
    control_value,
    seed=42,
    pseudocount=1.0,
    require_cell_cycle=True,
):
    """Run the standard validity scorecard for a perturbation-vs-control comparison.

    Returns a Scorecard with per-check PASS/FAIL/ARTIFACT verdicts plus an overall
    label. Use this BEFORE interpreting any module-level effect, as in paper4.5+ flow.

    overall_verdict semantics:
      PASS         all primary checks PASS (random_200 PASS, housekeeping PASS,
                   and if require_cell_cycle: cell_cycle PASS)
      ARTIFACT     random_200 FAIL → baseline shift confounds the analysis;
                   strongly suggests retrying with score_method='scanpy_corrected'
      INCONCLUSIVE perturbation didn't reach biology (cell_cycle FAIL); module-level
                   nulls may be uninformative at all
      MIXED        some checks PASS, some FAIL — interpret with caution
    """
    panels = _load_panels()
    thr = panels["thresholds"]

    rand_genes = random_200_panel(adata, seed=seed)
    rand_lfc = log2fc_perturbation_vs_control(
        adata, rand_genes, condition_col, perturbation_value, control_value, pseudocount,
    )
    hk_lfc = log2fc_perturbation_vs_control(
        adata, housekeeping_panel(), condition_col, perturbation_value, control_value, pseudocount,
    )
    cc_lfc = log2fc_perturbation_vs_control(
        adata, cell_cycle_panel(), condition_col, perturbation_value, control_value, pseudocount,
    )

    rand_v = _classify(rand_lfc, thr["random_200_max_abs_log2fc"], "lt")
    hk_v = _classify(hk_lfc, thr["housekeeping_max_abs_log2fc"], "lt")
    cc_v = _classify(cc_lfc, thr["cell_cycle_min_abs_log2fc"], "gt")

    checks = [
        CheckResult(
            name="random_200", value=rand_lfc,
            threshold=thr["random_200_max_abs_log2fc"], direction="lt",
            verdict=rand_v,
            notes=("baseline-shift artifact detector; FAIL ⇒ retry with "
                   "score_method='scanpy_corrected' (paper4.5 v1.2 fix)"),
        ),
        CheckResult(
            name="housekeeping", value=hk_lfc,
            threshold=thr["housekeeping_max_abs_log2fc"], direction="lt",
            verdict=hk_v,
            notes="technical normalization sanity",
        ),
        CheckResult(
            name="cell_cycle", value=cc_lfc,
            threshold=thr["cell_cycle_min_abs_log2fc"], direction="gt",
            verdict=cc_v,
            notes="positive control; FAIL ⇒ perturbation may not reach biology",
        ),
    ]

    if rand_v == "FAIL":
        overall = "ARTIFACT"
    elif require_cell_cycle and cc_v == "FAIL":
        overall = "INCONCLUSIVE"
    elif all(c.verdict == "PASS" for c in checks):
        overall = "PASS"
    else:
        overall = "MIXED"

    return Scorecard(overall_verdict=overall, checks=checks)


def scorecard(scorecard_obj_or_dict):
    """Pretty-print a Scorecard."""
    if hasattr(scorecard_obj_or_dict, "to_dict"):
        d = scorecard_obj_or_dict.to_dict()
    else:
        d = scorecard_obj_or_dict
    lines = [f"Overall: {d['overall_verdict']}", "Checks:"]
    for c in d["checks"]:
        sym = {"PASS": "✓", "FAIL": "✗", "ARTIFACT": "!", "INCONCLUSIVE": "?", "NA": "—"}.get(
            c["verdict"], "?"
        )
        op = "<" if c["direction"] == "lt" else ">"
        lines.append(
            f"  {sym} {c['name']:14s} value={c['value']:+.4f}  ({op} {c['threshold']:.2f})  → {c['verdict']}"
        )
        if c["notes"]:
            lines.append(f"      ↳ {c['notes']}")
    return "\n".join(lines)
