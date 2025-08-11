"""
A class representing a binding pocket in a protein structure.

Attributes:
    file_path (Optional[Path]): Path to the PDB file containing the pocket.
    block_type (str): Format of the block content ('pdb').
    block_content (str): String containing the pocket data.
"""

from dataclasses import dataclass, field
import io
import os
from pathlib import Path
import shutil
from typing import Any, Dict, List, Optional

from biotite.structure import AtomArrayStack
from biotite.structure.io.pdb import PDBFile
from deeporigin_molstar import ProteinViewer
import numpy as np
import pandas as pd
from tabulate import tabulate
from termcolor import colored

from deeporigin.drug_discovery.utilities.visualize import jupyter_visualization


@dataclass
class Pocket:
    """A class representing a binding pocket in a protein structure."""

    file_path: Optional[Path] = None
    block_type: str = ""
    block_content: str = ""
    color: str = "red"
    name: Optional[str] = None
    pdb_id: Optional[str] = None
    index: Optional[int] = 0
    props: Optional[Dict[str, Any]] = field(default_factory=dict)
    structure: Optional[np.ndarray] = None
    coordinates: Optional[np.ndarray] = None

    def __post_init__(self):
        if self.file_path is not None:
            self.load_structure(self.file_path)

    def load_structure(self, structure_file_path: str | Path) -> None:
        """
        Load a PDB structure from a file path into the `structure` attribute.

        Args:
            structure_file_path (str | Path): Path to the PDB file.
        """
        structure_file = PDBFile.read(str(structure_file_path))
        structure = structure_file.get_structure()

        if isinstance(structure, AtomArrayStack):
            self.structure = structure[0]
        else:
            self.structure = structure

        if self.structure is not None:
            self.coordinates = self.structure.coord

    def __repr__(self):
        # Basic info table
        basic_info = [
            ["Name", self.name],
            ["File", self.file_path],
            ["Color", self.color],
        ]

        # Properties table if available
        if self.props:
            properties = [
                ["Volume", f"{self.props.get('volume', 'N/A')} Å³"],
                ["Total SASA", f"{self.props.get('total_SASA', 'N/A')} Å²"],
                ["Polar SASA", f"{self.props.get('polar_SASA', 'N/A')} Å²"],
                [
                    "Polar/Apolar SASA ratio",
                    f"{self.props.get('polar_apolar_SASA_ratio', 'N/A')}",
                ],
                ["Hydrophobicity", f"{self.props.get('hydrophobicity', 'N/A')}"],
                ["Polarity", f"{self.props.get('polarity', 'N/A')}"],
                [
                    "Drugability score",
                    f"{self.props.get('drugability_score', 'N/A')}",
                ],
            ]
            return f"Pocket:\n{tabulate(basic_info, tablefmt='rounded_grid')}\n\nProperties:\n{tabulate(properties, tablefmt='rounded_grid')}"

        return f"Pocket:\n{tabulate(basic_info, tablefmt='rounded_grid')}"

    def get_center(self) -> np.ndarray:
        """
        Get the center of the pocket based on its coordinates.

        Returns:
            np.ndarray: A numpy array containing the center of the pocket.
        """
        if self.structure is None:
            raise ValueError("No structure loaded for this pocket")
        return self.structure.coord.mean(axis=0)

    @classmethod
    def from_pocket_finder_results(
        cls,
        pocket_finder_results_dir: str | Path,
    ) -> List["Pocket"]:
        """Create a list of Pocket objects from pocket finder results directory.

        Args:
            pocket_finder_results_dir: Directory containing pocket finder results
                with PDB files for each pocket and a CSV properties file.

        Returns:
            List of Pocket objects with properties from the CSV file.
        """
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
                file_path=pdb_file,
                name=pdb_file.stem,
                pdb_id=None,  # Will be set from CSV if available
                props={},  # Initialize empty properties dictionary
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
                pdb_filename = os.path.basename(pocket.file_path)

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

                            pocket.props[column] = value
        except Exception:
            # If there's an error reading the CSV, just return the pockets with basic info
            pass

        return pockets

    @classmethod
    def from_block(
        cls,
        block_content: str,
        block_type: str = "pdb",
        **kwargs: Any,  # noqa: ANN401
    ) -> "Pocket":
        """Create a Pocket instance from block content.

        Args:
            block_content (str): The content of the pocket structure.
            block_type (str): The format of the block content (default: "pdb").
            **kwargs: Additional arguments to pass to the Pocket constructor.

        Returns:
            Pocket: A new Pocket instance.


        """
        if not block_type:
            raise ValueError(
                "block_type must be provided when initializing with block_content."
            )

        pocket = cls(block_content=block_content, block_type=block_type, **kwargs)
        pocket._initialize_from_block()
        return pocket

    @classmethod
    def from_file(cls, file_path: str, **kwargs: Any) -> "Pocket":
        """
        Create a Pocket instance from a file.

        Args:
            file_path (str): Path to the pocket structure file.
            **kwargs: Additional arguments to pass to the Pocket constructor.

        Returns:
            Pocket: A new Pocket instance.


        """
        file_path = Path(file_path).absolute()
        if not file_path.exists():
            raise FileNotFoundError(f"The file {file_path} does not exist.")

        block_type = file_path.suffix.lstrip(".").lower()
        block_content = file_path.read_text()

        pocket = cls(
            file_path=file_path,
            block_type=block_type,
            block_content=block_content,
            **kwargs,
        )
        pocket._initialize_from_file()
        return pocket

    @classmethod
    def from_name(
        cls,
        name: str,
        **kwargs: Any,
    ) -> "Pocket":
        """
        Create a Pocket instance by searching for a file with the given name in the pockets directory.

        Args:
            name (str): Name of the pocket to search for.
            **kwargs: Additional arguments to pass to the Pocket constructor.

        Returns:
            Pocket: A new Pocket instance.


        """
        pockets_dir = Path(cls.get_directory())
        matching_files = list(pockets_dir.glob(f"{name}.*"))

        if not matching_files:
            raise FileNotFoundError(
                f"No pocket file found with name '{name}' in {pockets_dir}"
            )
        if len(matching_files) > 1:
            raise ValueError(
                f"Multiple files found with name '{name}' in {pockets_dir}"
            )

        return cls.from_file(matching_files[0], **kwargs)

    @classmethod
    def from_structure(
        cls,
        structure: np.ndarray,
        name: Optional[str] = None,
        **kwargs: Any,
    ) -> "Pocket":
        """
        Create a Pocket instance directly from a structure array.

        Args:
            structure (np.ndarray): The structure array.
            name (Optional[str]): Name for the pocket.
            **kwargs: Additional arguments to pass to the Pocket constructor.

        Returns:
            Pocket: A new Pocket instance.

        """
        pocket = cls(structure=structure, name=name, **kwargs)
        pocket._initialize_from_structure()
        return pocket

    def _initialize_from_block(self):
        """Initialize the pocket from block content."""
        try:
            self.structure = self.load_structure_from_block(
                self.block_content,
                self.block_type,
            )
            self._post_structure_initialization()

            # Save to file
            directory = Path(self.get_directory())
            self.file_path = directory / f"{self.name}.{self.block_type}"
            self.write_to_file(self.file_path)
        except Exception as e:
            raise RuntimeError(
                f"Failed to initialize pocket from block: {str(e)}"
            ) from e

    def _initialize_from_file(self):
        """Initialize the pocket from a file."""
        try:
            pocket_file_dir = self.get_directory()
            if str(pocket_file_dir) != str(self.file_path.parent):
                destination = Path(pocket_file_dir) / self.file_path.name
                shutil.copy2(self.file_path, destination)
                self.file_path = destination

            self.structure = self.load_structure_from_block(
                self.block_content, self.block_type
            )
            self._post_structure_initialization()
        except Exception as e:
            raise RuntimeError(
                f"Failed to initialize pocket from file: {str(e)}"
            ) from e

    def _initialize_from_structure(self):
        """Initialize the pocket from a structure array."""
        try:
            if isinstance(self.structure, AtomArrayStack):
                self.structure = self.structure[0]
            self._post_structure_initialization()

            # Save to file
            directory = Path(self.get_directory())
            self.file_path = directory / f"{self.name}.pdb"
            self.write_to_file(self.file_path)
        except Exception as e:
            raise RuntimeError(
                f"Failed to initialize pocket from structure: {str(e)}"
            ) from e

    def _post_structure_initialization(self):
        """Common initialization steps after structure is loaded."""
        if self.structure is None:
            raise ValueError("Structure could not be loaded.")

        if self.name is None:
            if self.file_path:
                self.name = self.file_path.stem
            else:
                self.name = "Unknown_Pocket"
                directory = Path(self.get_directory())
                num = len(list(directory.glob(f"{self.name}*")))
                self.name = f"{self.name}_{num + 1}"

        self.coordinates = self.structure.coord

    def load_structure_from_block(self, block_content: str, block_type: str):
        """
        Load a pocket structure from block content.

        Parameters:
        - block_content (str): String containing the pocket data.
        - block_type (str): Format of the block content ('pdb').

        Returns:
        - AtomArray: Loaded structure.

        Raises:
        - ValueError: If the block type is unsupported.
        """
        if block_type == "pdb":
            pdb_file = PDBFile.read(io.StringIO(block_content))
            structure = pdb_file.get_structure()
        else:
            raise ValueError(f"Unsupported block type: {block_type}")
        return structure

    def write_to_file(self, output_path: str, output_format: str = "pdb"):
        """
        Write the current state of the structure to a PDB file.

        Parameters:
        - file_path (str): Path where the pocket structure will be written.

        """

        def write_to_pdb_file(structure, output_path):
            pdb_file = PDBFile()
            pdb_file.set_structure(structure)
            pdb_file.write(output_path)

        try:
            path = Path(output_path)
            if not path.parent.exists():
                path.parent.mkdir(parents=True, exist_ok=True)

            if path.suffix.lower() != ".pdb":
                raise NotImplementedError("convert_file not implemented yet")
                # with tempfile.NamedTemporaryFile(delete=True) as temp:
                #     write_to_pdb_file(self.structure, temp.name)
                #     convert_file("pdb", temp.name, output_format, output_path)
            else:
                write_to_pdb_file(self.structure, output_path)

        except Exception as e:
            raise RuntimeError(
                f"Failed to write structure to file {output_path}: {str(e)}"
            ) from e

    @jupyter_visualization
    def show(self):
        """show the pocket in a jupyter notebook"""

        pocket_paths = [str(self.file_path)]
        pocket_names = ["Name: " + self.name + " | " + self.pocket_props()]

        viewer = ProteinViewer("", format="pdb")

        pocket_config = viewer.get_pocket_visualization_config()
        # Ensure self.index is within the bounds of surface_colors
        if self.index >= len(pocket_config.surface_colors):
            self.index = 0  # Default to the first color if out of bounds

        pocket_config.surface_colors = [pocket_config.surface_colors[self.index]]
        print(
            "\n|\n".join(
                colored("■", pocket_config.surface_colors[0]) + " " + pocket_name
                for pocket_name in pocket_names
            )
        )

        return viewer.render_protein_with_pockets(
            pocket_paths=pocket_paths, pocket_config=pocket_config
        )

    def pocket_props(self):
        """get the properties of the pocket"""

        properties_line = ""
        if self.props:
            properties_line = (
                f"Volume: {self.props.get('volume', 'N/A')}Å³ | "
                f"Drugability score: {self.props.get('drugability_score', 'N/A')}"
            )
        return properties_line

    def __str__(self):
        properties_line = ""
        if self.props:
            properties_line = (
                f"  Volume: {self.props.get('volume', 'N/A')}Å³, "
                f"Total SASA: {self.props.get('total_SASA', 'N/A')}, "
                f"Polar SASA: {self.props.get('polar_SASA', 'N/A')}, "
                f"Polar/Apolar SASA ratio: {self.props.get('polar_apolar_SASA_ratio', 'N/A')}, "
                f"Hydrophobicity: {self.props.get('hydrophobicity', 'N/A')}, "
                f"Polarity: {self.props.get('polarity', 'N/A')}, "
                f"Drugability score: {self.props.get('drugability_score', 'N/A')}"
            )

        return (
            f"Pocket:\n  Name: {self.name}\n{properties_line}  Block type: {self.block_type}\n"
            "Available Fields: {block_type, block_content, file_path, name, coordinates}"
        )

    @staticmethod
    def get_directory() -> str:
        """
        Generates and ensures the existence of a directory for pockets.

        Returns:
            str: The path to the pockets directory.
        """
        home_dir = Path.home()
        pockets_base_dir = home_dir / ".deeporigin" / "pockets"
        pockets_base_dir.mkdir(parents=True, exist_ok=True)

        return str(pockets_base_dir)

    def update_coordinates(self, coords: np.ndarray):
        """update coordinates of the pocket structure"""

        self.structure.coord = coords
        self.coordinates = coords
