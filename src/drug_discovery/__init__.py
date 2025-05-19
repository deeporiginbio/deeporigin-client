"""
Drug Discovery Module

This module provides tools and utilities for drug discovery workflows, including
molecule manipulation, protein-ligand interactions, and computational chemistry
calculations.
"""

# example data
import asyncio  # noqa: F401
from importlib.resources import path

import nest_asyncio  # we need this for asyncio  # noqa: F401

from . import chemistry

__all__ = ["chemistry", "Complex", "Protein", "Ligand", "Pocket"]

with path("deeporigin.data.brd", "brd.pdb") as file_path:
    EXAMPLE_DATA_DIR = file_path.parent


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
    elif name == "Pocket":
        from .structures.pocket import Pocket

        return Pocket
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")


def __dir__():
    return __all__ + ["EXAMPLE_DATA_DIR"]
