"""Exact and conditional audits for fixed MERCI-style score inputs."""
from .analytic import analytic
from .conditioned import evaluate, solve
from .inputs import validate_manifest
from .sampler import GateSampler

__all__ = ["GateSampler", "analytic", "evaluate", "solve", "validate_manifest"]
__version__ = "0.1.1"
