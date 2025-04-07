import io
import shutil
import tempfile

import numpy as np
from pathlib import Path
from typing import Optional, List
from termcolor import colored

from biotite.structure.io.pdb import PDBFile
from biotite.structure import AtomArrayStack, centroid

from deeporigin_molstar import ProteinViewer

from deeporigin.src.constants import WORKING_DIR
from deeporigin.src.utilities.logging import DEFAULT_LOGGER
from deeporigin.src.utilities.visualize import jupyter_visualization
from deeporigin.src.utilities.conversions import convert_file


class Pocket:
    def __init__(
        self,
        file_path: str = "",
        block_type: str = "",
        block_content: str = "",
        color="red",
        name=None,
        index: Optional[int] = 0,
        props: Optional[dict] = None,
    ):
        self.color = color
        self.index = index
        self.props = props

        self.name = name
        self.file_path = None
        self.structure = None
        self.coordinates = None

        # Determine block_type from file extension if not provided
        file_path_obj = Path(file_path).absolute() if file_path else None
        extension = file_path_obj.suffix.lower() if file_path_obj else ""
        if not block_type and extension:
            block_type = extension.lstrip(".")  # Remove the leading dot
        self.block_type = block_type.lower()
        self.block_content = block_content

        # Ensure only one source is provided
        sources_provided = sum(bool(x) for x in [file_path, block_content])
        if sources_provided != 1:
            raise ValueError(
                "Please provide exactly one of file_path or block_content."
            )

        from_block = False
        try:
            if file_path:
                self.file_path = Path(file_path).absolute()
                if not self.file_path.exists():
                    raise FileNotFoundError(
                        f"The file {self.file_path} does not exist."
                    )

                if not self.block_type:
                    self.block_type = self.file_path.suffix.lstrip(".").lower()

                self.block_content = self.file_path.read_text()

                pocket_file_dir = self.get_directory()
                if str(pocket_file_dir) != str(self.file_path.parent):
                    try:
                        destination = Path(pocket_file_dir) / self.file_path.name
                        shutil.copy2(self.file_path, destination)
                        self.file_path = destination
                    except Exception as e:
                        DEFAULT_LOGGER.log_error(
                            f"Failed to copy file to destination: {str(e)}"
                        )
                        raise
            elif block_content:
                self.block_content = block_content
                if not self.block_type:
                    raise ValueError(
                        "block_type must be provided when initializing with block_content."
                    )
                from_block = True
                pocket_file_dir = self.get_directory()

            if self.block_content:
                if self.block_type not in ["pdb"]:
                    raise ValueError(
                        f"Only pdb file formats are supported (given {self.block_type})"
                    )
                self.structure = self.load_structure_from_block(
                    self.block_content, self.block_type
                )

            if self.structure is None:
                raise ValueError("Structure could not be loaded.")

            # Type checking for AtomArrayStack
            if isinstance(self.structure, AtomArrayStack):
                self.structure = self.structure[0]

            DEFAULT_LOGGER.log_info(
                f"Loaded structure from {self.file_path if self.file_path else 'block content'}. Selected structure index: {0}"
            )

            if self.name is None:
                if self.file_path:
                    self.name = self.file_path.stem
                else:
                    self.name = "Unknown_Pocket"
                    directory = Path(pocket_file_dir)
                    num = len(list(directory.glob(f"{self.name}*")))
                    self.name = f"{self.name}_{num + 1}"

            self.coordinates = self.structure.coord

            if from_block:
                directory = Path(pocket_file_dir)
                self.file_path = directory / f"{self.name}.{self.block_type}"
                self.write_to_file(self.file_path)

        except Exception as e:
            DEFAULT_LOGGER.log_error(f"Failed to initialize pocket: {str(e)}")
            raise

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

    @staticmethod
    def load_structure(structure_file_path: str):
        """
        Load a PDB structure from a file path.

        Parameters:
        - structure_file_path (str): Path to the PDB file.

        Returns:
        - AtomArray: Loaded structure.
        """
        structure_file = PDBFile.read(structure_file_path)
        structure = structure_file.get_structure()
        return structure

    def write_to_file(self, output_path: str, output_format: str = "pdb"):
        """
        Write the current state of the structure to a PDB file.

        Parameters:
        - file_path (str): Path where the pocket structure will be written.

        Example:
        ```python
        pocket.write_to_file('/path/to/output.pdb')
        ```
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
                with tempfile.NamedTemporaryFile(delete=True) as temp:
                    write_to_pdb_file(self.structure, temp.name)
                    convert_file("pdb", temp.name, output_format, output_path)
            else:
                write_to_pdb_file(self.structure, output_path)
            DEFAULT_LOGGER.log_info(f"Current structure written to {output_path}.")

        except Exception as e:
            DEFAULT_LOGGER.log_error(
                f"Failed to write structure to file {output_path}: {str(e)}"
            )

    @jupyter_visualization
    def visualize(self):
        pocket_paths = [str(self.file_path)]
        pocket_names = ["Name: " + self.name + " | " + self.pocket_props()]

        viewer = ProteinViewer("", format="pdb")

        pocket_config = viewer.get_pocket_visualization_config()
        # Ensure self.index is within the bounds of surface_colors
        if self.index >= len(pocket_config.surface_colors):
            # Log a warning or adjust index appropriately
            DEFAULT_LOGGER.log_warning(
                f"Index {self.index} is out of bounds for surface_colors. Resetting to 0."
            )
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
        properties_line = ""
        if self.props:
            properties_line = (
                f"Volume: {self.props.get('volume', 'N/A')}Å³ | "
                f"Drugability score: {self.props.get('drugability_score', 'N/A')}"
            )
        return properties_line

    def __repr__(self):
        properties_line = ""
        if self.props:
            properties_line = (
                f"  Volume: {self.props.get('volume', 'N/A')}Å³, "
                f"Total SASA: {self.props.get('total_SASA', 'N/A')} "
                f"Polar SASA: {self.props.get('polar_SASA', 'N/A')} "
                f"Polar/Apolar SASA ratio: {self.props.get('polar_apolar_SASA_ratio', 'N/A')} "
                f"Hydrophobicity: {self.props.get('hydrophobicity', 'N/A')} "
                f"Polarity: {self.props.get('polarity', 'N/A')} "
                f"Drugability score: {self.props.get('drugability_score', 'N/A')}"
            )

        return (
            f"Pocket:\n  Name: {self.name}\n{properties_line}  Block type: {self.block_type}\n"
            "Available Fields: {block_type, block_content, file_path, name, coordinates}"
        )

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

    def get_center(self) -> Optional[List[float]]:
        """
        Get the center of the ligand based on its coordinates.

        Returns:
        - list: The center coordinates of the ligand.
        - None: If coordinates are not available.

        Example:
        ```python
        center = ligand.get_center()
        print(center)  # Output: [1.23, 4.56, 7.89]
        ```
        """
        if self.coordinates is None:
            DEFAULT_LOGGER.log_warning("Coordinates are not available for this Pocket.")
            return None
        center = self.coordinates.mean(axis=0)
        DEFAULT_LOGGER.log_info(f"Calculated center coordinates: {center.tolist()}")
        return [float(x) for x in center.tolist()]

    @staticmethod
    def get_directory() -> str:
        """
        Generates and ensures the existence of a directory for pockets.

        Returns:
            str: The path to the pockets directory.
        """
        pockets_base_dir = Path(WORKING_DIR) / "pockets"
        pockets_base_dir.mkdir(parents=True, exist_ok=True)

        return str(pockets_base_dir)


    def update_coordinates(self, coords: np.ndarray):
        self.structure.coord = coords
        self.coordinates = coords
        DEFAULT_LOGGER.log_info("Pocket coordinates has been inplaced updated.")