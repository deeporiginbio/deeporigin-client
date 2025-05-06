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
from .structures import Ligand, Pocket, Protein

__all__ = ["chemistry", "Complex", "Protein", "Ligand", "Pocket"]


# Lazy imports for better performance
def __getattr__(name):
    if name == "Complex":
        from deeporigin.drug_discovery.complex import Complex

        return Complex

    raise AttributeError(f"module {__name__} has no attribute {name}")


with path("deeporigin.data.brd", "brd.pdb") as file_path:
    EXAMPLE_DATA_DIR = file_path.parent
