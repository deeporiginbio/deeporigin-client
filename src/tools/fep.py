"""module to work FEP calculations"""

from dataclasses import dataclass, fields
from pathlib import Path
from typing import Optional, Union

from deeporigin import chemistry


@dataclass
class Ligand:
    """class to represent a ligand (typically backed by a SDF file)"""

    file: Union[str, Path]
    smiles_string: Optional[str] = None

    def __post_init__(self):
        """generates a SMILES if it doesn't exist"""

        if self.smiles_string is None:
            self.smiles_string = chemistry.sdf_to_smiles(self.file)

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
