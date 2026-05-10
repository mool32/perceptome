"""Module co-variation network + per-module heterogeneity analyses."""

from .network import compute_network, compare_networks
from .heterogeneity import module_heterogeneity

__all__ = ["compute_network", "compare_networks", "module_heterogeneity"]
