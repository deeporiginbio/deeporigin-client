import tempfile
from copy import deepcopy
from pathlib import Path
from typing import List, Optional

import numpy as np
from deeporigin_molstar import DockingViewer as DockingMolstarViewer

from deeporigin.src.utilities.alignments import StructureAligner, create_bounding_box
from deeporigin.src.utilities.logging import DEFAULT_LOGGER
from deeporigin.src.utilities.utils import (
    calculate_box_dimensions,
    calculate_box_min_max,
)
from deeporigin.src.utilities.visualize import jupyter_visualization

from .ligand import Ligand
from .pocket import Pocket
from .protein import Protein


class PocketData:
    def __init__(
        self,
        pocket: List | Pocket | Ligand,
        xref_protein: Protein,
        box_size: Optional[List] = None,
        fit_box: Optional[bool] = None,
        padding: Optional[int] = 2,
    ):
        if box_size and fit_box:
            raise ValueError("Cannot specify both box_size and [fit_box and padding]")

        if not box_size and not fit_box:
            DEFAULT_LOGGER.log_warning(
                "No box size or fit box specified. Defaulting to 24Å box."
            )
            box_size = [24, 24, 24]

        self.__pocket = deepcopy(pocket)
        if isinstance(self.__pocket, list):
            if len(self.__pocket) != 3:
                raise ValueError(
                    "If pocket is a list, it must have 3 elements (x, y, z)"
                )

            for i in self.__pocket:
                if not isinstance(
                    i, (int, float, np.int32, np.int64, np.float32, np.float64)
                ):
                    raise ValueError(
                        "All elements in pocket must be integers or floats"
                    )

            self.__aligner = None

            if fit_box:
                raise ValueError("Cannot fit box to list of coordinates")

        elif isinstance(pocket, (Pocket, Ligand)):
            self.__aligner = StructureAligner()
            self.__aligner.calculate_pca(self.__pocket.coordinates)

            pocket_transformed_coords = self.transform(self.__pocket.coordinates)
            self.__pocket.update_coordinates(pocket_transformed_coords)
        else:
            raise ValueError("Invalid pocket type")

        self.__padding = padding
        self.__fit_box = fit_box
        self.__xref_protein = xref_protein
        self.__box_center = self._get_center(self.__pocket)
        self.__box_min_coords, self.__box_max_coords, self.__box_size = (
            self.calculate_box_params(box_size)
        )

    @classmethod
    def create_from_residues(
        cls,
        xref_protein: Protein,
        residue_ids: List[str],
        box_size: Optional[List] = None,
        fit_box: Optional[bool] = None,
        padding: Optional[int] = 2,
    ):
        center, warning, _ = xref_protein.get_center_by_residues(residue_ids)
        if warning:
            DEFAULT_LOGGER.log_warning(warning)

        return cls(
            pocket=center,
            xref_protein=xref_protein,
            box_size=box_size,
            fit_box=fit_box,
            padding=padding,
        )

    @property
    def pocket(self):
        return self.__pocket

    @property
    def xref_protein(self):
        return self.__xref_protein

    @property
    def padding(self):
        return self.__padding

    @property
    def fit_box(self):
        return self.__fit_box

    @property
    def box_size(self):
        return self.__box_size

    @property
    def box_center(self):
        return self.__box_center

    @property
    def aligner(self):
        return self.__aligner

    @property
    def box_min_coords(self):
        return self.__box_min_coords

    @property
    def box_max_coords(self):
        return self.__box_max_coords

    def _get_center(self, pocket: str | Pocket | List[float]):
        """
        Convert pocket and box_size to lists of floats.

        Args:
            pocket: Pocket info.
            box_size: Box dimensions.

        Returns:
            (pocket_center, box_size) as float lists.
        """
        if isinstance(pocket, list):
            pocket = [float(x) for x in pocket]
        elif isinstance(pocket, (Pocket, Ligand)):
            pocket = pocket.get_center()
        else:
            raise ValueError("Invalid pocket type.")

        return pocket

    def transform(self, coords: np.ndarray) -> np.ndarray:
        """
        Transform coordinates to the PCA-aligned space.

        Args:
            coords: Coordinates to transform.

        Returns:
            Transformed coordinates.
        """
        if self.aligner is None:
            return coords

        return self.aligner.align_structure(coords)

    def inverse_transform(self, coords: np.ndarray) -> np.ndarray:
        """
        Inverse transform coordinates from the PCA-aligned space.

        Args:
            coords: Coordinates to inverse transform.

        Returns:
            Inverse transformed coordinates.
        """
        if self.aligner is None:
            return coords

        return self.aligner.restore_structure(coords)

    def match_protein(self, protein: Protein) -> bool:
        """
        Check if the protein matches the xref protein.

        Args:
            protein: Protein to check.

        Returns:
            True if the proteins match, False otherwise.
        """
        if protein != self.xref_protein:
            raise ValueError("Provided protein does not match xref protein.")

    def calculate_box_params(self, box_size: Optional[List] = None):
        if self.fit_box:
            result = create_bounding_box(
                self.pocket, padding=self.padding, around_ligand=True
            )
            box_min_coords, box_max_coords = (
                list(result["min_coords"]),
                list(result["max_coords"]),
            )
        else:
            box_min_coords, box_max_coords = calculate_box_min_max(
                box_center=self.box_center, box_dimensions=box_size
            )

        if self.fit_box:
            box_size = calculate_box_dimensions(box_min_coords, box_max_coords)

        return box_min_coords, box_max_coords, box_size

    def from_xyz(self):
        if isinstance(self.pocket, (Ligand, Pocket)):
            return False
        elif isinstance(self.pocket, list):
            return True
        else:
            raise ValueError("Invalid pocket type")

    @jupyter_visualization
    def show_box(
        self,
        protein: Protein = None,
        raise_for_protein_mismatch: bool = True,
    ) -> str:
        """
        Display protein-ligand complex or bounding box visualization.

        Args:
            protein: Protein object.
            structure: Ligand, Pocket or coordinates.
            fit_box: Fit box around structure.
            padding: Padding if fitting.
            box_size: Box dimensions if not fitting.

        Returns:
            HTML representation of the visualization.
        """
        if protein is None:
            protein = self.xref_protein

        if raise_for_protein_mismatch:
            self.match_protein(protein)

        protein = deepcopy(protein)
        aligned_protein_coords = self.transform(protein.structure.coord)
        protein.update_coordinates(aligned_protein_coords)

        if not self.from_xyz():
            with tempfile.TemporaryDirectory() as temp_dir:
                tmp_protein_file = Path(temp_dir) / "protein.pdb"
                tmp_structure_file = Path(temp_dir) / "structure.sdf"

                protein.write_to_file(str(tmp_protein_file))
                self.pocket.write_to_file(str(tmp_structure_file), output_format="sdf")

                html = DockingMolstarViewer().render_ligand_with_bounding_box(
                    protein_data=str(tmp_protein_file),
                    protein_format="pdb",
                    ligand_data=str(tmp_structure_file),
                    ligand_format="sdf",
                    box={"min": self.box_min_coords, "max": self.box_max_coords},
                )
        else:
            protein_file_path = str(protein.file_path)
            protein_format = getattr(protein, "block_type", "pdb")
            html = DockingMolstarViewer().render_bounding_box(
                protein_data=protein_file_path,
                protein_format=protein_format,
                box_center=self.box_center,
                box_size=self.box_size,
            )

        return html
