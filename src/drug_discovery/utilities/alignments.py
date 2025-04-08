from typing import Optional

import biotite.structure as struc
import numpy as np
from sklearn.decomposition import PCA

from deeporigin.src.utilities.logging import DEFAULT_LOGGER
from deeporigin.src.utilities.utils import calculate_box_min_max


def calculate_bounding_box(structure_coord, padding=0.0):
    """
    Calculate the bounding box dimensions and center for a given structure.

    Parameters:
    - structure (AtomArray): The structure containing atom coordinates.
    - padding (float): Padding to add around the bounding box.

    Returns:
    - tuple: min_coords, max_coords, dimensions, and center.
    """
    min_coords = np.min(structure_coord, axis=0) - padding
    max_coords = np.max(structure_coord, axis=0) + padding
    dimensions = max_coords - min_coords
    center = 0.5 * (max_coords + min_coords)

    return min_coords, max_coords, dimensions, center


def calculate_fixed_bounding_box(structure_coord, box_size=24.0):
    """
    Calculate a fixed bounding box centered on the structure.

    Parameters:
    - structure (AtomArray): The structure containing atom coordinates.
    - box_size (float): The size of the bounding box in angstroms (default is 24.0).

    Returns:
    - tuple: min_coords, max_coords, dimensions, and center.
    """
    # Calculate the centroid of the structure
    center = np.mean(structure_coord, axis=0)

    # Calculate half the box size
    half_size = box_size / 2.0

    # Define min and max coordinates based on the center and half_size
    min_coords = center - half_size
    max_coords = center + half_size

    # Dimensions of the box
    dimensions = max_coords - min_coords

    return min_coords, max_coords, dimensions, center


def create_bounding_box_atoms(min_coords, max_coords, dimensions):
    """
    Create atoms at the corners of the bounding box.

    Parameters:
    - min_coords (np.ndarray): Minimum coordinates of the bounding box.
    - max_coords (np.ndarray): Maximum coordinates of the bounding box.
    - dimensions (np.ndarray): Dimensions of the bounding box.

    Returns:
    - AtomArray: An array of atoms at the bounding box corners.
    """
    atoms = []
    atom_index = 1

    # Create an atom at the minimum corner
    atom = struc.Atom(
        min_coords,
        chain_id="A",
        res_id=1,
        res_name="BOX",
        atom_name=f"XE{atom_index}",
        element="XE",
    )
    atoms.append(atom)
    atom_index += 1

    # Create atoms at the corners formed by adding dimensions to min_coords
    for i in range(3):
        coord = min_coords.copy()
        coord[i] += dimensions[i]
        atom = struc.Atom(
            coord,
            chain_id="A",
            res_id=1,
            res_name="BOX",
            atom_name=f"XE{atom_index}",
            element="XE",
        )
        atoms.append(atom)
        atom_index += 1

    # Create an atom at the maximum corner
    atom = struc.Atom(
        max_coords,
        chain_id="A",
        res_id=1,
        res_name="BOX",
        atom_name=f"XE{atom_index}",
        element="XE",
    )
    atoms.append(atom)
    atom_index += 1

    # Create atoms at the corners formed by subtracting dimensions from max_coords
    for i in range(3):
        coord = max_coords.copy()
        coord[i] -= dimensions[i]
        atom = struc.Atom(
            coord,
            chain_id="A",
            res_id=1,
            res_name="BOX",
            atom_name=f"XE{atom_index}",
            element="XE",
        )
        atoms.append(atom)
        atom_index += 1

    # Convert the list of atoms to an AtomArray
    atom_array = struc.array(atoms)
    return atom_array


def create_bounding_box(
    ligand, padding=0.0, output_file=None, around_ligand=False, box_size=20.0
):
    """
    Calculates the bounding box and optionally creates a PDB file with atoms at the corners.

    Parameters:
    - ligand (Ligand): Ligand structure to create a bounding box for.
    - padding (float): Padding to add to the bounding box dimensions.
    - output_file (str): If provided, saves the atoms to this PDB file.

    Returns:
    - dict: Contains min_coords, max_coords, dimensions, center, and optionally 'atom_array' if output_file is provided.
    """
    structure_coord = ligand.coordinates
    if around_ligand:
        min_coords, max_coords, dimensions, center = calculate_bounding_box(
            structure_coord, padding
        )
    else:
        min_coords, max_coords, dimensions, center = calculate_fixed_bounding_box(
            structure_coord, box_size
        )

    result = {
        "min_coords": min_coords,
        "max_coords": max_coords,
        "dimensions": dimensions,
        "center": center,
    }

    if output_file:
        atom_array = create_bounding_box_atoms(min_coords, max_coords, dimensions)
        struc.io.save_structure(output_file, atom_array)
        result["atom_array"] = atom_array
        print(f"Bounding box atoms saved to {output_file}")

    return result


def save_bounding_box(center, box_size, output_file):
    """
    Save atoms at the corners of the bounding box to a PDB file.

    Parameters:
    - center (np.ndarray): Box center.
    - box_size (np.ndarray): Dimensions of the bounding box.
    - output_file (str): File to save the atoms to.
    """
    min_coords, max_coords = calculate_box_min_max(center, box_size)
    atom_array = create_bounding_box_atoms(min_coords, max_coords, box_size)
    struc.io.save_structure(output_file, atom_array)


class StructureAligner:
    """
    A class to handle structure alignment using PCA transformations.
    Maintains PCA state and provides methods for alignment and restoration.
    """

    def __init__(self):
        """Initialize the StructureAligner with no initial PCA."""
        self.pca: Optional[PCA] = None
        self._components_fixed: bool = False

    def calculate_pca(self, coords: np.ndarray) -> None:
        """
        Calculate and store PCA components using provided coordinates.

        Parameters:
        - coords: np.ndarray containing coordinates for PCA calculation

        Raises:
        - ValueError: If coordinates are None
        - Exception: For general PCA calculation failures
        """
        try:
            if coords is None:
                raise ValueError("Coordinates are None.")

            self.pca = PCA(n_components=3)
            self.pca.fit(coords)

            components = self.pca.components_
            if np.dot(np.cross(components[0], components[1]), components[2]) < 0:
                self.pca.components_ *= np.array([1, 1, -1])

            self._components_fixed = True
            DEFAULT_LOGGER.log_info("PCA components calculated and stored.")

        except Exception as e:
            DEFAULT_LOGGER.log_error(f"PCA calculation failed: {str(e)}")
            raise

    @property
    def is_fitted(self) -> bool:
        """Check if PCA has been calculated."""
        return self.pca is not None and self._components_fixed

    def align_structure(self, coords: np.ndarray) -> np.ndarray:
        """
        Align coordinates using the stored PCA components.

        Parameters:
        - coords: np.ndarray to align

        Returns:
        - np.ndarray: Aligned coordinates

        Raises:
        - ValueError: If PCA hasn't been calculated or coordinates are None
        - Exception: For general alignment failures
        """
        if not self.is_fitted:
            raise ValueError(
                "PCA components haven't been calculated. Call calculate_pca first."
            )

        try:
            if coords is None:
                raise ValueError("Coordinates are None.")

            # Transform coordinates
            aligned_coords = self.pca.transform(coords)
            DEFAULT_LOGGER.log_info("Coordinates aligned using PCA.")

            return aligned_coords

        except Exception as e:
            DEFAULT_LOGGER.log_error(f"Alignment failed: {str(e)}")
            raise

    def restore_structure(self, coords: np.ndarray) -> np.ndarray:
        """
        Restore coordinates from PCA-aligned space back to original space.

        Parameters:
        - coords: np.ndarray to restore

        Returns:
        - np.ndarray: Restored coordinates in original space

        Raises:
        - ValueError: If PCA hasn't been calculated or coordinates are None
        - Exception: For general restoration failures
        """
        if not self.is_fitted:
            raise ValueError(
                "PCA components haven't been calculated. Call calculate_pca first."
            )

        try:
            if coords is None:
                raise ValueError("Coordinates are None.")

            # Inverse transform coordinates
            restored_coords = self.pca.inverse_transform(coords)
            DEFAULT_LOGGER.log_info("Coordinates restored from PCA space.")

            return restored_coords

        except Exception as e:
            DEFAULT_LOGGER.log_error(f"Restoration failed: {str(e)}")
            raise
