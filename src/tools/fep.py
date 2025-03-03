"""module to work FEP calculations"""

import importlib.resources
import json
import os
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Optional, Union

import pandas as pd
from beartype import beartype
from box import Box
from deeporigin import chemistry
from deeporigin.tools.toolkit import _ensure_columns, _ensure_database
from IPython.display import HTML, display


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

    # this ID keeps track of whether it is uploaded to deep origin or not
    _do_id: Optional[str] = None

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

    # this ID keeps track of whether it is uploaded to deep origin or not
    _do_id: Optional[str] = None

    def __post_init__(self):
        self.file = Path(self.file)
        if self.name is None:
            self.name = self.file.name


@beartype
def _ensure_dbs_for_fep():
    """ensure that there are databases for FEP on the data hub

    Ensures the following DBs:

    - ligand (list of ligand files)
    - protein (list of protein files)
    - ABFE (ligand, protein, step, job_id, output, results, delta_g)
    - RBFE (ligand1, ligand2, protein, step, job_id, output, results, delta_g)

    """

    # ligands
    database = _ensure_database("Ligands")
    required_columns = [
        dict(name="Ligand", type="file"),
    ]
    _ensure_columns(
        database=database,
        required_columns=required_columns,
    )

    # proteins
    database = _ensure_database("Proteins")
    required_columns = [
        dict(name="Protein", type="file"),
    ]
    _ensure_columns(
        database=database,
        required_columns=required_columns,
    )

    # ABFE
    database = _ensure_database("ABFE")
    required_columns = [
        dict(name="Protein", type="text"),
        dict(name="Ligand", type="text"),
        dict(name="Step", type="text"),
        dict(name="JobID", type="text"),
        dict(name="OutputFile", type="file"),
        dict(name="ResultFile", type="file"),
        dict(name="delta_g", type="file"),
    ]
    _ensure_columns(
        database=database,
        required_columns=required_columns,
    )


@dataclass
class FEP:
    """class for FEP simulations"""

    ligands: list[Ligand]
    protein: Protein

    @classmethod
    def from_dir(cls, directory: str) -> "FEP":
        """initialize an FEP class given some files in a directory

        Args:
            directory (str): directory containing ligand and protein files"""

        sdf_files = sorted(
            [
                os.path.join(directory, f)
                for f in os.listdir(directory)
                if f.lower().endswith(".sdf")
            ]
        )
        ligands = [Ligand(sdf_file) for sdf_file in sdf_files]

        pdb_files = [
            os.path.join(directory, f)
            for f in os.listdir(directory)
            if f.lower().endswith(".pdb")
        ]

        if len(pdb_files) != 1:
            raise ValueError(
                f"Expected exactly one PDB file in the directory, but found {len(pdb_files)}."
            )
        protein_file = pdb_files[0]
        protein = Protein(protein_file)

        # Create the ABFE instance
        fep = cls(
            ligands=ligands,
            protein=protein,
        )

        return fep

    def show_ligands(self):
        """show ligands in a dataframe"""

        # convert SMILES to aligned images
        smiles_list = [ligand.smiles_string for ligand in self.ligands]
        id_list = [ligand._do_id for ligand in self.ligands]
        file_list = [os.path.basename(ligand.file) for ligand in self.ligands]

        images = chemistry.smiles_list_to_base64_png_list(smiles_list)

        data = {
            "ID": id_list,
            "File": file_list,
            "Ligand": images,
        }

        df = pd.DataFrame(data)

        # Use escape=False to allow the <img> tags to render as images
        display(HTML(df.to_html(escape=False)))

    def connect(self):
        """connects the local instantiation of the simulation to Deep Origin"""

        # ensure that ligands are uploaded
        _ensure_dbs_for_fep()
