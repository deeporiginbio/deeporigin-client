"""
Drug Discovery Module

This module provides tools and utilities for drug discovery workflows, including
molecule manipulation, protein-ligand interactions, and computational chemistry
calculations.
"""

# example data
import asyncio  # noqa: F401
from importlib.resources import files

import nest_asyncio  # we need this for asyncio  # noqa: F401

from . import chemistry

__all__ = ["chemistry", "Complex", "Protein", "Ligand", "Pocket", "LigandSet"]

DATA_DIR = files("deeporigin.data")
BRD_DATA_DIR = DATA_DIR / "brd"


def __getattr__(name):
    if name == "Complex":
        from .complex import Complex

        return Complex
    elif name == "Protein":
        from .structures.protein import Protein

        return Protein
    elif name == "Ligand":
        from .structures.ligand import Ligand

        return Ligand
    elif name == "LigandSet":
        from .structures.ligand import LigandSet

        return LigandSet
    elif name == "Pocket":
        from .structures.pocket import Pocket

        return Pocket
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return __all__ + ["BRD_DATA_DIR", "DATA_DIR"]
