"""module to work FEP calculations"""

import importlib.resources
import json
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Optional, Union

from beartype import beartype
from box import Box
from deeporigin import chemistry


class PrettyDict(Box):
    """A dict subclass with a custom pretty-print representation."""

    def __repr__(self):
        """pretty print a dict"""
        return json.dumps(
            dict(self),
            indent=2,
            ensure_ascii=False,
        )

    def _repr_html_(self):
        """pretty print a dict"""
        self.__repr__()


@beartype
def _load_params(step: str) -> Box:
    """load default values for abfe end to end run"""

    with importlib.resources.open_text("deeporigin.json", f"{step}.json") as f:
        return PrettyDict(json.load(f))


@dataclass
class Ligand:
    """class to represent a ligand (typically backed by a SDF file)"""

    file: Union[str, Path]
    smiles_string: Optional[str] = None
    n_molecules: Optional[int] = None

    def __post_init__(self):
        """generates a SMILES if it doesn't exist"""

        # check that there's only one molecule here
        if self.n_molecules is None:
            if chemistry.count_molecules_in_sdf_file(self.file) > 1:
                raise ValueError(
                    "Too many molecules. Expected a single molecule in the SDF file, but got multiple"
                )
            self.n_molecules = 1

        if self.smiles_string is None:
            smiles_string = chemistry.sdf_to_smiles(self.file)
            if len(smiles_string) > 1:
                raise ValueError("Expected a single SMILES strings, but got multiple")
            self.smiles_string = smiles_string[0]

    def _repr_pretty_(self, p, cycle):
        """pretty print a ligand"""

        if cycle:
            p.text("Ligand(...)")
        else:
            p.text("Ligand(")

            with p.group(2, "\n  ", "\n"):
                all_fields = fields(self)
                for idx, field in enumerate(all_fields):
                    value = getattr(self, field.name)
                    p.text(f"{field.name}: {value!r}")
                    # Only add a breakable if this isn't the last field.
                    if idx < len(all_fields) - 1:
                        p.breakable()
            p.text(")")


@dataclass
class Protein:
    """class to represent a protein (typically backed by a PDB file)"""

    file: Union[str, Path]
    name: Optional[str] = None

    def __post_init__(self):
        self.file = Path(self.file)
        if self.name is None:
            self.name = self.file.name
