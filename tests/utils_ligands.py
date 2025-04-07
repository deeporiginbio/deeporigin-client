"""Shared test fixtures and setup code for drug discovery tests"""

import os

base_path = os.path.join(os.path.dirname(__file__), "fixtures")

ligands = [
    {
        "file": os.path.join(base_path, "42-ligands.sdf"),
        "n_ligands": 42,
        "name_by_property": "Compound",
    },
    {
        "file": os.path.join(base_path, "ligands-brd-all.sdf"),
        "n_ligands": 8,
        "name_by_property": "_Name",
    },
    {
        "file": os.path.join(base_path, "brd-7.sdf"),
        "n_ligands": 1,
        "name_by_property": "_Name",
    },
]


bad_ligands = [
    {
        "file": "missing.sdf",
        "smiles_string": None,
    },
    {
        "file": None,
        "smiles_string": None,
    },
]
