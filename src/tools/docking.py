"""this module contains tools for docking on Deep Origin"""

import os
from dataclasses import dataclass
from typing import Optional

import pandas as pd
from beartype import beartype
from deeporigin import chemistry
from deeporigin.data_hub import api
from deeporigin.tools.toolkit import _ensure_columns, _ensure_database
from deeporigin.utils.core import hash_strings


@dataclass
class Docking:
    """class to run docking on Deep Origin"""

    protein: chemistry.Protein
    smiles_strings: list[str]

    _proteins_db: Optional[dict] = None
    _docking_db: Optional[dict] = None

    _all_smiles_hash: Optional[str] = None

    def _repr_pretty_(self, p, cycle):
        """pretty print a Docking object"""

        if cycle:
            p.text("Docking(...)")
        else:
            p.text("Docking(")

            p.text(f"protein={self.protein.name}")
            p.text(f" with {len(self.smiles_strings)} ligands")
            p.text(")")

    def show_results(self):
        """show results of bulk Docking run in a table, rendering 2D structures of molecules"""

        df = self.get_results()

        from IPython.display import HTML, display

        smiles_list = list(df["SMILES"])
        images = chemistry.smiles_list_to_base64_png_list(smiles_list)
        df["Structure"] = images
        df.drop("SMILES", axis=1, inplace=True)
        display(HTML(df.to_html(escape=False)))
