"""module that contains classes for various chemical structures, including Ligands, Proteins, Complexes, and Pockets"""

import os
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Optional, Tuple

import numpy as np
import pandas as pd
from beartype import beartype
from biotite.structure import AtomArrayStack
from biotite.structure.io.pdb import PDBFile
from tabulate import tabulate

from deeporigin.drug_discovery import chemistry
from deeporigin.exceptions import DeepOriginException
from deeporigin.functions.pocket_finder import find_pockets


@dataclass
class Ligand:
    """Class to represent a ligand (typically backed by a SDF file)"""

    file: Optional[str | Path] = None
    smiles_string: Optional[str] = None

    # this ID keeps track of whether it is uploaded to deep origin or not
    _do_id: Optional[str] = None

    # this stores user-defined properties
    properties: Optional[dict] = None

    def __post_init__(self):
        """post init tasks"""

        if self.file is not None and not os.path.exists(self.file):
            raise DeepOriginException(f"File {self.file} does not exist")

        # we require either a SMILES string or a file
        if self.file is None and self.smiles_string is None:
            raise DeepOriginException("Must specify either a file or a SMILES string")

        # read user-defined properties
        if self.file is not None:
            self.properties = chemistry.read_sdf_properties(self.file)

            # check that there's only one molecule here
            if chemistry.count_molecules_in_sdf_file(self.file) > 1:
                raise ValueError(
                    "Too many molecules. Expected a single molecule in the SDF file, but got multiple"
                )

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

    def show(self):
        """show a ligand in a Jupyter notebook using molstar"""

        if self.file is not None:
            # backed by SDF file. use a 3D viewer

            chemistry.show_molecules_in_sdf_file(self.file)
        else:
            # only SMILES. use a 2D viewer
            img = chemistry.smiles_to_base64_png(self.smiles_string)

            from IPython.display import HTML, display

            display(HTML(img))

    @classmethod
    def from_smiles(cls, smiles: str) -> "Ligand":
        """create a ligand from a SMILES string"""
        return cls(smiles_string=smiles)

    @classmethod
    def from_csv(
        cls,
        *,
        file: str | Path,
        smiles_column: str,
        properties_columns: list[str] = None,
    ) -> list["Ligand"]:
        """create a list of ligands from a CSV file

        Args:
            file: Path to CSV file
            smiles_column: Column name containing SMILES strings
            properties_columns: List of column names to extract as properties

        Returns:
            List of Ligand objects
        """

        # Read the CSV file
        df = pd.read_csv(file)

        # Validate column existence
        if smiles_column not in df.columns:
            raise ValueError(f"SMILES column '{smiles_column}' not found in CSV file")

        # Create empty list to store ligands
        ligands = []

        # Process each row
        for _, row in df.iterrows():
            smiles = row[smiles_column]

            # Skip empty SMILES
            if pd.isna(smiles) or not smiles.strip():
                continue

            # Extract properties if columns were specified
            properties = None
            if properties_columns:
                properties = {}
                for col in properties_columns:
                    if col in df.columns:
                        properties[col] = row[col]
                    else:
                        # Skip non-existent columns with a warning
                        print(f"Warning: Property column '{col}' not found in CSV file")

            # Create ligand and add to list
            ligand = cls(smiles_string=smiles, properties=properties)
            ligands.append(ligand)

        return ligands


@dataclass
class Pocket:
    """Class to represent a pocket in a protein"""

    file: Optional[str | Path] = None
    name: Optional[str] = None
    pdb_id: Optional[str] = None
    properties: Optional[dict] = None
    color: str = "red"
    structure = None

    def __post_init__(self):
        if self.file is not None:
            self.load_structure(self.file)

    def load_structure(self, structure_file_path: str):
        """
        Load a PDB structure from a file path.

        Args:
            - structure_file_path (str): Path to the PDB file.

        Returns:
            AtomArray: Loaded structure.
        """
        structure_file = PDBFile.read(structure_file_path)

        structure = structure_file.get_structure()

        if isinstance(structure, AtomArrayStack):
            self.structure = structure[0]

    def __repr__(self):
        # Basic info table
        basic_info = [
            ["Name", self.name],
            ["File", self.file],
            ["Color", self.color],
        ]

        # Properties table if available
        if self.properties:
            properties = [
                ["Volume", f"{self.properties.get('volume', 'N/A')} Å³"],
                ["Total SASA", f"{self.properties.get('total_SASA', 'N/A')} Å²"],
                ["Polar SASA", f"{self.properties.get('polar_SASA', 'N/A')} Å²"],
                [
                    "Polar/Apolar SASA ratio",
                    f"{self.properties.get('polar_apolar_SASA_ratio', 'N/A')}",
                ],
                ["Hydrophobicity", f"{self.properties.get('hydrophobicity', 'N/A')}"],
                ["Polarity", f"{self.properties.get('polarity', 'N/A')}"],
                [
                    "Drugability score",
                    f"{self.properties.get('drugability_score', 'N/A')}",
                ],
            ]
            return f"Pocket:\n{tabulate(basic_info, tablefmt='rounded_grid')}\n\nProperties:\n{tabulate(properties, tablefmt='rounded_grid')}"

        return f"Pocket:\n{tabulate(basic_info, tablefmt='rounded_grid')}"

    def get_center(self) -> Optional[list[float]]:
        """
        Get the center of the pocket based on its coordinates.

        Returns:
            - list of floats: The center coordinates of the pocket.

        """

        center = self.structure.coord.mean(axis=0)

        return center

        # return [float(x) for x in center.tolist()]

    @classmethod
    def from_pocket_finder_results(
        cls,
        pocket_finder_results_dir: str | Path,
    ) -> "list[Pocket]":
        """Create a list of Pocket objects from pocket finder results directory.

        Args:
            pocket_finder_results_dir: Directory containing pocket finder results
                with PDB files for each pocket and a CSV properties file.

        Returns:
            List of Pocket objects with properties from the CSV file.
        """

        from pathlib import Path

        import pandas as pd

        # Convert to Path object for consistent handling
        results_dir = Path(pocket_finder_results_dir)

        # Find all PDB files in the directory
        pdb_files = list(results_dir.glob("*.pdb"))

        # Colors to cycle through for pockets
        colors = [
            "red",
            "green",
            "blue",
            "yellow",
            "orange",
            "gray",
            "purple",
            "cyan",
            "magenta",
            "lime",
        ]

        # Create Pocket objects from each PDB file
        pockets = []
        for idx, pdb_file in enumerate(pdb_files):
            # Create a Pocket object with the file path and name (without extension)
            pocket = cls(
                file=str(pdb_file),
                name=pdb_file.stem,
                pdb_id=None,  # Will be set from CSV if available
                properties={},  # Initialize empty properties dictionary
                color=colors[idx % len(colors)],  # Cycle through colors
            )
            pockets.append(pocket)

        # Find the CSV properties file
        csv_files = list(results_dir.glob("*.csv"))
        if not csv_files:
            # If no CSV file, return pockets with just file and name
            return pockets

        # Use the first CSV file found
        properties_file = csv_files[0]

        try:
            # Read CSV file using pandas
            df = pd.read_csv(properties_file)

            # Map CSV columns to Pocket properties
            # Using 'pocket_file' column to match with PDB file names
            for pocket in pockets:
                # Get the filename without extension to match with CSV
                pdb_filename = os.path.basename(pocket.file)

                # Try to find a row in the CSV that matches this pocket file
                pocket_row = df[df["pocket_file"] == pdb_filename]

                if not pocket_row.empty:
                    # Add all properties from the CSV to the pocket's properties dictionary
                    for column in df.columns:
                        if column != "pocket_file":  # Skip the file column
                            value = pocket_row[column].iloc[0]

                            # Convert NumPy types to Python primitive types
                            if isinstance(value, np.integer):
                                value = int(value)
                            elif isinstance(value, np.floating):
                                value = float(value)
                            elif isinstance(value, np.bool_):
                                value = bool(value)
                            elif pd.isna(value):
                                value = None

                            pocket.properties[column] = value
        except Exception as e:
            # Log the error but still return the pockets with basic info
            import logging

            logging.warning(f"Error reading properties from CSV: {e}")

        return pockets


@dataclass
class Protein:
    """Class to represent a protein (typically backed by a PDB file)"""

    file: Optional[str | Path] = None
    name: Optional[str] = None
    pdb_id: Optional[str] = None

    # this ID keeps track of whether it is uploaded to deep origin or not
    _do_id: Optional[str] = None

    def __post_init__(self):
        if self.pdb_id is not None:
            self.file = chemistry.download_protein(self.pdb_id)

        self.file = Path(self.file)
        if self.name is None:
            self.name = self.file.name

    @beartype
    def show(
        self,
        pockets: Optional[list[Pocket]] = None,
        sdf_file=None,
    ):
        """visualize the protein in a Jupyter notebook using molstar"""

        from deeporigin_molstar import JupyterViewer, ProteinViewer

        if pockets is None and sdf_file is None:
            protein_viewer = ProteinViewer(
                data=str(self.file),
                format="pdb",
            )
            html_content = protein_viewer.render_protein()

            JupyterViewer.visualize(html_content)
        elif pockets is not None and sdf_file is None:
            pocket_surface_alpha: float = 0.7
            protein_surface_alpha: float = 0.1

            protein_viewer = ProteinViewer(data=str(self.file), format="pdb")
            pocket_paths = [pocket.file for pocket in pockets]

            # Retrieve and customize pocket visualization configuration
            pocket_config = protein_viewer.get_pocket_visualization_config()
            pocket_config.surface_alpha = pocket_surface_alpha

            protein_config = protein_viewer.get_protein_visualization_config()
            protein_config.style_type = "cartoon"
            protein_config.surface_alpha = protein_surface_alpha
            pocket_config.surface_colors = [pocket.color for pocket in pockets]

            # Render the protein with pockets
            html_content = protein_viewer.render_protein_with_pockets(
                pocket_paths=pocket_paths,
                pocket_config=pocket_config,
                protein_config=protein_config,
            )
            from deeporigin_molstar import JupyterViewer

            JupyterViewer.visualize(html_content)

        elif sdf_file is not None:
            from deeporigin_molstar import DockingViewer, JupyterViewer

            docking_viewer = DockingViewer()
            html_content = docking_viewer.render_with_seperate_crystal(
                protein_data=str(self.file),
                protein_format="pdb",
                ligands_data=[sdf_file],
                ligand_format="sdf",
            )

            JupyterViewer.visualize(html_content)

    @beartype
    def find_pockets(
        self,
        pocket_count: int = 5,
        pocket_min_size: int = 30,
    ) -> list[Pocket]:
        """find pockets in the protein"""
        results_dir = find_pockets(
            self.file,
            pocket_count=pocket_count,
            pocket_min_size=pocket_min_size,
        )

        return Pocket.from_pocket_finder_results(results_dir)

    def dock(
        self,
        *,
        ligand: Ligand,
        pocket: Pocket,
    ):
        from deeporigin.functions import docking

        docked_ligand_sdf_file = docking.dock(
            protein=self,
            ligand=ligand,
            pocket=pocket,
        )

        return docked_ligand_sdf_file


@beartype
def ligands_to_dataframe(ligands: list[Ligand]):
    """convert a list of ligands to a pandas dataframe"""

    import pandas as pd

    smiles_list = [ligand.smiles_string for ligand in ligands]
    id_list = [ligand._do_id for ligand in ligands]
    file_list = [
        os.path.basename(ligand.file) if ligand.file is not None else None
        for ligand in ligands
    ]

    data = {
        "Ligand": smiles_list,
        "ID": id_list,
        "File": file_list,
    }

    # find all the common properties in all ligands
    common_keys = set.intersection(
        *(set(ligand.properties.keys()) for ligand in ligands)
    )
    for key in common_keys:
        data[key] = [ligand.properties[key] for ligand in ligands]

    df = pd.DataFrame(data)

    return df


@beartype
def show_ligands(ligands: list[Ligand]):
    """show ligands in the FEP object in a dataframe. This function visualizes the ligands using core-aligned 2D visualizations.

    Args:
        ligands (list[Ligand]): list of ligands

    """

    df = ligands_to_dataframe(ligands)

    # convert SMILES to aligned images
    images = chemistry.smiles_list_to_base64_png_list(df["Ligand"].tolist())
    df["Ligand"] = images

    # Use escape=False to allow the <img> tags to render as images
    from IPython.display import HTML, display

    display(HTML(df.to_html(escape=False)))
