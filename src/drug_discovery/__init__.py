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

from deeporigin.drug_discovery.complex import Complex

from . import chemistry
from .structures import Ligand, Pocket, Protein

__all__ = ["chemistry", "Complex", "Protein", "Ligand", "Pocket"]


with path("deeporigin.data.brd", "brd.pdb") as file_path:
    EXAMPLE_DATA_DIR = file_path.parent
