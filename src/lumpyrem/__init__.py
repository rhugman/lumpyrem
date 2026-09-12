"""lumpyrem -- lumped-parameter recharge modelling.

A Python port of the LUMPREM/LUMPREM2 Fortran models, adding array-based
forcing and output, a modern API, and compiled-speed execution.

Phase 0 (reference oracle and scaffolding) is complete; the model kernel
arrives in Phase 1.  See docs/conversion-plan.md.
"""

__version__ = "0.0.1"
__all__ = ["__version__"]
