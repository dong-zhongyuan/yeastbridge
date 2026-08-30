"""Target-agnostic CRC candidate-universe and real-perturbation baseline scan."""

from .core import ScanError, build_target_universe, pareto_fronts, run_scan

__all__ = ["ScanError", "build_target_universe", "pareto_fronts", "run_scan"]
