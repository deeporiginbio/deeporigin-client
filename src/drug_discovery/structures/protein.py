"""
Protein Module

This module encapsulates the Protein class, which is responsible for managing and manipulating
protein structures in computational biology workflows. It provides functionalities to load protein data
from various sources, preprocess structures, handle ligands, and visualize protein-ligand interactions.

Key Features:
- Protein Class: Manages protein data loaded from PDB IDs, file paths, or direct block content.
  Handles structure selection, filtering of HETATM records, and extraction of metal and cofactor residues.
- Initialization: Allows initializing a Protein instance using a PDB ID, local file, or block content,
  ensuring only one source is provided and properly loading the structure.
- Preparation Methods: Includes methods like `prepare`, `remove_hetatm`, `remove_resnames`, and
  `remove_water` to preprocess and clean protein structures for downstream analysis.
- Ligand Management: Provides methods to select and manage ligands within the protein structure,
  enabling extraction and manipulation of ligand residues.
- Visualization: Integrates with the ProteinViewer to render interactive visualizations of protein
  structures within Jupyter Notebooks, facilitating intuitive analysis and presentation of protein-ligand complexes.
- Utility Functions: Offers static and internal helper methods for downloading protein data, loading structures,
  and creating new Protein instances with modified structures.

Dependencies:
- Biotite for structure handling
- RDKit for cheminformatics
- Plotly for interactive visualizations
- External tools for enhanced functionality

Usage Example:
```python
# Initialize with a PDB ID
protein = Protein.from_pdb_id("1ABC")

# Initialize with a file
protein = Protein.from_file("/path/to/protein.pdb")

# Select a specific chain
chain_a = protein.select_chain('A')

# Remove water molecules
protein_no_water = protein.remove_water()

# Visualize the protein structure
protein.show()
```
"""

from collections import defaultdict
from dataclasses import dataclass, field
import io
import os
from pathlib import Path
import tempfile
from typing import Optional, Tuple

from beartype import beartype
from biotite.database.rcsb import fetch
from biotite.structure import filter_solvent
from biotite.structure.geometry import centroid
from biotite.structure.io.pdb import PDBFile
from deeporigin_molstar import DockingViewer, JupyterViewer, ProteinViewer
import numpy as np

from deeporigin.drug_discovery.constants import METALS, STATE_DUMP_PATH
from deeporigin.drug_discovery.external_tools.utils import (
    generate_html_output,
    get_protein_info_dict,
)
from deeporigin.functions.pocket_finder import find_pockets

from .entity import Entity
from .ligand import Ligand
from .pocket import Pocket


@dataclass
class Protein(Entity):
    """A class representing a protein structure with various manipulation and analysis capabilities."""

    # Core attributes
    structure: np.ndarray = field(repr=False)
    name: str
    file_path: Optional[Path] = None
    pdb_id: Optional[str] = None
    info: Optional[dict] = None
    atom_types: Optional[np.ndarray] = None
    block_type: str = "pdb"
    block_content: Optional[str] = None

    _remote_path_base = "entities/proteins/"

    @classmethod
    def from_pdb_id(cls, pdb_id: str, struct_ind: int = 0) -> "Protein":
        """
        Create a Protein instance from a PDB ID.

        Args:
            pdb_id (str): PDB ID of the protein to download.
            struct_ind (int): Index of the structure to select if multiple are present.

        Returns:
            Protein: A new Protein instance.

        Raises:
            ValueError: If the PDB ID is invalid or the structure cannot be loaded.
            RuntimeError: If the download fails.
        """
        try:
            file_path = Path(cls.download_protein_by_pdb_id(pdb_id)).absolute()
            block_content = file_path.read_text()
            structure = cls.load_structure_from_block(block_content, "pdb")
            structure = cls.select_structure(structure, struct_ind)

            return cls(
                structure=structure,
                name=pdb_id,
                file_path=file_path,
                pdb_id=pdb_id,
                info=get_protein_info_dict(pdb_id),
                atom_types=structure.atom_name,
                block_content=block_content,
            )
        except Exception as e:
            raise RuntimeError(
                f"Failed to create Protein from PDB ID {pdb_id}: {str(e)}"
            ) from e

    @classmethod
    def from_file(cls, file_path: str, struct_ind: int = 0) -> "Protein":
        """
        Create a Protein instance from a file.

        Args:
            file_path (str): Path to the protein PDB file.
            struct_ind (int): Index of the structure to select if multiple are present.

        Returns:
            Protein: A new Protein instance.

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If the structure cannot be loaded.
            RuntimeError: If the file cannot be read or processed.
        """
        try:
            file_path = Path(file_path).absolute()
            if not file_path.exists():
                raise FileNotFoundError(f"The file {file_path} does not exist.")

            block_type = file_path.suffix.lstrip(".").lower()
            block_content = file_path.read_text()
            structure = cls.load_structure_from_block(block_content, block_type)
            structure = cls.select_structure(structure, struct_ind)

            return cls(
                structure=structure,
                name=file_path.stem,
                file_path=file_path,
                atom_types=structure.atom_name,
                block_type=block_type,
                block_content=block_content,
            )
        except FileNotFoundError:
            raise
        except Exception as e:
            raise RuntimeError(
                f"Failed to create Protein from file {file_path}: {str(e)}"
            ) from e

    @staticmethod
    def load_structure_from_block(block_content: str, block_type: str) -> np.ndarray:
        """Load a protein structure from block content."""
        if block_type in ["pdb", "pdbqt"]:
            pdb_file = PDBFile.read(io.StringIO(block_content))
            structure = pdb_file.get_structure()
        else:
            raise ValueError(f"Unsupported block type: {block_type}")
        return structure

    @staticmethod
    def select_structure(structure: np.ndarray, index: int) -> np.ndarray:
        """Select a specific structure by index."""
        if index < 0 or index >= len(structure):
            raise ValueError(
                f"Invalid structure index {index}. Total structures: {len(structure)}"
            )
        return structure[index]

    @staticmethod
    def download_protein_by_pdb_id(pdb_id: str, save_dir: str = "") -> str:
        """Download a PDB structure by its PDB ID from RCSB."""
        if save_dir == "":
            save_dir = Protein.get_directory()

        pdb_id = pdb_id.lower()
        save_dir_path = Path(save_dir)
        save_dir_path.mkdir(parents=True, exist_ok=True)

        file_path = save_dir_path / f"{pdb_id}.pdb"
        if not file_path.exists():
            try:
                fetch(pdb_id, "pdb", save_dir_path)
            except Exception as e:
                raise RuntimeError(f"Failed to download PDB {pdb_id}: {str(e)}") from e

        return str(file_path)

    @staticmethod
    def get_directory() -> str:
        """Get the directory for storing protein files."""
        home_dir = Path.home()
        proteins_dir = home_dir / ".deeporigin" / "proteins"
        proteins_dir.mkdir(parents=True, exist_ok=True)
        return str(proteins_dir)

    @beartype
    def dock(
        self,
        *,
        ligand: Ligand,
        pocket: Pocket,
    ) -> str:
        """Dock a ligand into a specific pocket of the protein.

        This method performs molecular docking of a ligand into a specified pocket
        of the protein structure. It uses the Deep Origin docking to
        generate a 3D structure of the docked ligand.

        Args:
            ligand (Ligand): The ligand to dock into the protein pocket.
            pocket (Pocket): The specific pocket in the protein where the ligand
                should be docked.

        Returns:
            str: Path to the SDF file containing the docked ligand structure.
        """
        from deeporigin.functions import docking

        docked_ligand_sdf_file = docking.dock(
            protein=self,
            ligand=ligand,
            pocket=pocket,
        )

        return docked_ligand_sdf_file

    @property
    def coordinates(self):
        return self.structure.coord

    # @beartype
    # def prepare(self, model_loops: bool = False, pdb_id: str = "") -> "Protein":
    #     """
    #     Prepares the protein by calling the 'prepare' function with the specified protein path, PDB ID, and extension.
    #     It extracts metal and cofactor residue names from the protein structure and passes them to the 'prepare' function.

    #     Returns:
    #         protein (Protein): The prepared Protein object.
    #     Raises:
    #         Exception: If the preparation of the protein fails.
    #     """
    #     pdb_id = pdb_id if pdb_id else self.pdb_id
    #     if model_loops and not pdb_id:
    #         raise ValueError("PDB ID must be provided to model loops.")

    #     metal_resnames, cofactor_resnames = self.extract_metals_and_cofactors()
    #     metals_to_keep = [
    #         resname for resname in metal_resnames if resname.upper() in METALS
    #     ]

    #     response = prepare(
    #         protein_path=self.file_path,
    #         protein_pdb_id=pdb_id,
    #         protein_extension=self.block_type,
    #         metal_resnames=metals_to_keep,
    #         cofactor_resnames=cofactor_resnames,
    #         model_loops=model_loops,
    #     )
    #     if not response["prepared_protein_content"]:
    #         raise Exception("Failed to prepare protein.")

    #     protein_dir = Path(self.file_path).parent
    #     base_name = (
    #         Path(self.file_path).stem if self.file_path else "modified_structure"
    #     )
    #     new_file_name = protein_dir / f"{base_name}_prep.pdb"

    #     intermediate_protein = Protein(
    #         block_content=response["prepared_protein_content"], block_type="pdb"
    #     )
    #     intermediate_protein.write_to_file(str(new_file_name))

    #     protein = Protein(file_path=new_file_name)
    #     protein.pdb_id = self.pdb_id

    #     return protein

    def _filter_hetatm_records(
        self, exclude_water: bool = True, keep_resnames: Optional[list[str]] = None
    ):
        """
        Internal method to filter HETATM records, optionally excluding water molecules and keeping specified residues.

        Parameters:
        - exclude_water (bool): Whether to exclude water molecules (default: True).
        - keep_resnames (Optional[list[str]]): List of residue names to keep (e.g., metal ions, cofactors).

        Returns:
        - AtomArray: Filtered HETATM records from the structure.
        """
        hetatm_records = self.structure[self.structure.hetero]
        res_names_upper = np.char.upper(hetatm_records.res_name)

        if exclude_water:
            water_residue_names = ["HOH", "WAT"]
            water_residue_names_upper = [name.upper() for name in water_residue_names]
            hetatm_records = hetatm_records[
                ~np.isin(res_names_upper, water_residue_names_upper)
            ]
            res_names_upper = np.char.upper(hetatm_records.res_name)

        if keep_resnames:
            keep_resnames_upper = [name.upper() for name in keep_resnames]
            hetatm_records = hetatm_records[
                np.isin(res_names_upper, keep_resnames_upper)
            ]

        return hetatm_records

    def _filter_chain_records(self, chain_ids: Optional[list[str]] = None):
        """
        Filter chain records based on chain IDs.

        Args:
        Parameters:
        - chain_ids (Optional[List[str]]): List of chain IDs to filter. If None or contains "ALL", all chains are returned.

        Returns:
        - AtomArray: Filtered chain records from the structure.
        """
        if chain_ids is None or "ALL" in chain_ids:
            return self.structure
        else:
            return self.structure[np.isin(self.structure.chain_id, chain_ids)]

    def list_chain_names(self) -> list[str]:
        """
        List all unique chain IDs in the protein structure.

        Returns:
            list[str]: A list of unique chain IDs.
        """
        chain_records = self._filter_chain_records()
        chain_ids = np.unique(chain_records.chain_id)
        return list(chain_ids)

    def list_hetero_names(self, exclude_water=True) -> list[str]:
        """
        List all unique hetero residue names in the protein structure.

        Args:
            exclude_water (bool): Whether to exclude water molecules from the list.

        Returns:
            list[str]: A list of unique ligand residue names (excluding water).
        """
        hetatm_records = self._filter_hetatm_records(exclude_water=exclude_water)
        ligand_res_names = np.unique(hetatm_records.res_name)
        return list(ligand_res_names)

    def select_chain(self, chain_id: str) -> Optional["Protein"]:
        """
        Select a specific chain by its ID and return a new Protein object.

        Parameters:
        - chain_id (str): Chain ID to select.

        Returns:
        - Protein: A new Protein object containing the selected chain.

        Raises:
        - ValueError: If the chain ID is not found.

        Example:
        ```python
        chain_a = protein.select_chain('A')
        ```
        """
        chain_records = self._filter_chain_records(chain_ids=[chain_id])
        if len(chain_records) > 0:
            return self._create_new_protein_with_structure(
                chain_records, suffix=f"_chain_{chain_id}"
            )
        else:
            raise ValueError(f"Chain {chain_id} not found.")

    def select_chains(self, chain_ids: list[str]) -> "Protein":
        """
        Select specific chains from the protein structure.

        Args:
            chain_ids (list[str]): List of chain IDs to select.
        """
        chain_records = self._filter_chain_records(chain_ids=chain_ids)
        if len(chain_records) == 0:
            raise ValueError(f"No chains found for the provided chain IDs: {chain_ids}")
        return self._create_new_protein_with_structure(
            chain_records, suffix=f"_chains_{'_'.join(chain_ids)}"
        )

    @beartype
    def find_pockets(
        self,
        pocket_count: int = 5,
        pocket_min_size: int = 30,
    ) -> list[Pocket]:
        """Find potential binding pockets in the protein structure.

        This method analyzes the protein structure to identify cavities or pockets
        that could potentially serve as binding sites for ligands. It uses the
        Deep Origin pocket finding algorithm to detect and characterize these pockets.

        Args:
            pocket_count (int, optional): Maximum number of pockets to identify.
                Defaults to 5.
            pocket_min_size (int, optional): Minimum size of pockets to consider,
                measured in cubic Angstroms. Defaults to 30.

        Returns:
            list[Pocket]: A list of Pocket objects, each representing a potential
                binding site in the protein. Each Pocket object contains:
                - The 3D structure of the pocket
                - Properties such as volume, surface area, hydrophobicity, etc.
                - Visualization parameters (color, etc.)

        Examples:
            >>> protein = Protein(file="protein.pdb")
            >>> pockets = protein.find_pockets(pocket_count=3, pocket_min_size=50)
            >>> for pocket in pockets:
            ...     print(f"Pocket: {pocket.name}, Volume: {pocket.properties.get('volume')} Å³")
        """
        results_dir = find_pockets(
            self.file_path,
            pocket_count=pocket_count,
            pocket_min_size=pocket_min_size,
        )

        return Pocket.from_pocket_finder_results(results_dir)

    @beartype
    def remove_hetatm(
        self,
        keep_resnames: Optional[list[str]] = None,
        remove_metals: Optional[list[str]] = None,
    ) -> None:
        """
        Remove HETATM records from the protein structure, with options to retain specified residues or exclude certain metals.

        Parameters
        ----------
        keep_resnames : Optional[list[str]]
            A list of residue names (strings) to keep in the structure even if they are HETATM records.
        exclude_metals : Optional[list[str]]
            A list of metal names (strings) to exclude from removal. These metals will be retained in the structure.

        Notes
        -----
        - By default, a predefined list of metals is considered for removal unless specified in `exclude_metals`.
        - If `keep_resnames` is provided, those residues (along with any metals not excluded) will be retained even if they are HETATM records.
        - The method updates the current protein object in place.

        Example:
        ```python
            protein = Protein(structure)
            protein.remove_hetatm(keep_resnames=['HOH'], exclude_metals=['ZN'])
        ```
        """
        metals = METALS
        if remove_metals:
            exclude_metals_upper = [metal.upper() for metal in remove_metals]
            metals = list(set(METALS) - set(exclude_metals_upper))

        if not metals and not keep_resnames:
            self.structure = self.structure[~self.structure.hetero]
        else:
            keep_resnames_upper = (
                [res.upper() for res in keep_resnames] if keep_resnames else []
            )
            keep_resnames_upper.extend(metals)
            keep_resnames_set = list(set(keep_resnames_upper))

            hetatm_to_keep = self._filter_hetatm_records(
                keep_resnames=keep_resnames_set
            )
            hetatm_indices_to_keep = np.isin(
                self.structure.res_id, hetatm_to_keep.res_id
            )
            self.structure = self.structure[
                ~self.structure.hetero | hetatm_indices_to_keep
            ]

    @beartype
    def remove_resnames(
        self,
        exclude_resnames: Optional[list[str]] = None,
    ) -> None:
        """
        Remove specific residue names from the protein structure in place.

        Args:
            exclude_resnames (Optional[list[str]]): List of residue names to exclude.
        """
        if exclude_resnames is not None:
            b_resn = np.isin(self.structure.res_name, exclude_resnames)
            self.structure = self.structure[~b_resn]

    def remove_water(self) -> None:
        """
        Remove water molecules from the protein structure in place.

        Example:
        ```python
        protein.remove_water()
        ```
        """
        self.structure = self.structure[~filter_solvent(self.structure)]

    def extract_metals_and_cofactors(self) -> Tuple[list[str], list[str]]:
        """
        Extract metal ions and cofactors from the protein structure.

        Returns:
            Tuple[list[str], list[str]]
        """
        hetatm_records = self.structure[self.structure.hetero]
        water_residue_names = ["HOH", "WAT"]
        hetatm_records = hetatm_records[
            ~np.isin(hetatm_records.res_name, water_residue_names)
        ]

        metal_elements = {
            "AC",
            "AG",
            "AL",
            "AM",
            "AS",
            "AU",
            "B",
            "BA",
            "BE",
            "BH",
            "BI",
            "BK",
            "CA",
            "CD",
            "CE",
            "CF",
            "CM",
            "CN",
            "CS",
            "CU",
            "DB",
            "DS",
            "DY",
            "ER",
            "ES",
            "EU",
            "FE",
            "FM",
            "FR",
            "GA",
            "GD",
            "GE",
            "HF",
            "HG",
            "HO",
            "HS",
            "K",
            "LA",
            "LI",
            "LR",
            "LU",
            "MD",
            "MG",
            "MN",
            "MO",
            "MT",
            "NA",
            "NB",
            "ND",
            "NI",
            "NO",
            "NP",
            "OS",
            "PA",
            "TA",
            "PM",
            "PO",
            "PR",
            "PT",
            "PU",
            "RA",
            "RB",
            "RE",
            "RF",
            "RG",
            "RH",
            "RU",
            "SB",
            "SC",
            "SG",
            "SI",
            "SM",
            "SN",
            "SR",
            "TB",
            "TC",
            "TE",
            "TH",
            "TI",
            "TL",
            "TM",
            "U",
            "V",
            "W",
            "YB",
            "ZN",
            "ZR",
            "CO",
            "CR",
            "IN",
            "IR",
            "PB",
            "PD",
        }

        residue_groups = defaultdict(list)
        for atom in hetatm_records:
            key = (atom.chain_id, atom.res_id, atom.ins_code)
            residue_groups[key].append(atom)

        metal_resnames = set()
        cofactor_resnames = set()
        for _, atoms in residue_groups.items():
            res_name = atoms[0].res_name.strip().upper()
            is_metal = all(
                atom.element.strip().upper() in metal_elements for atom in atoms
            )
            if is_metal:
                metal_resnames.add(res_name)
            else:
                cofactor_resnames.add(res_name)

        metal_resnames = list(metal_resnames)
        cofactor_resnames = list(cofactor_resnames)

        return metal_resnames, cofactor_resnames

    def _create_new_protein_with_structure(
        self, new_structure, suffix: str = "_modified"
    ) -> "Protein":
        """
        Helper method to create a new Protein object with a modified structure.
        Writes the modified structure to a new file and creates a new Protein object from that file.

        Parameters:
        - new_structure: The modified structure.
        - suffix (str): A suffix to append to the new file name (default: "_modified").

        Returns:
        - Protein: A new Protein object created from the newly written structure file.

        Raises:
        - Exception: If writing the new structure fails.
        """
        base_name = self.file_path.stem if self.file_path else "modified_structure"
        new_file_name = f"{base_name}{suffix}.pdb"
        parent_dir = (
            self.file_path.parent if self.file_path else Path(tempfile.gettempdir())
        )
        new_file_path = parent_dir / new_file_name

        if new_file_path.exists():
            os.remove(new_file_path)

        try:
            pdb_file = PDBFile()
            pdb_file.set_structure(new_structure)
            pdb_file.write(str(new_file_path))
            return Protein.from_file(str(new_file_path))
        except Exception as e:
            raise RuntimeError(
                f"Failed to create new Protein with modified structure: {str(e)}"
            ) from e

    def to_pdb(self, file_path: str):
        """
        Write the protein structure to a PDB file.

        Args:
            file_path (str): Path where the PDB file will be written.

        Example:
        ```python
        protein.to_pdb('/path/to/output.pdb')
        ```
        """
        try:
            pdb_file = PDBFile()
            pdb_file.set_structure(self.structure)
            pdb_file.write(file_path)
        except Exception as e:
            raise RuntimeError(
                f"Failed to write structure to file {file_path}: {str(e)}"
            ) from e

    @beartype
    def _dump_state(self) -> str:
        """Dump the current protein state to a fixed location in the user's home directory.

        Returns:
            str: Path to the state dump file containing the protein structure.
        """
        # Create the .deeporigin directory if it doesn't exist
        STATE_DUMP_PATH.parent.mkdir(exist_ok=True)

        # Use the constant file path
        self.to_pdb(str(STATE_DUMP_PATH))
        return str(STATE_DUMP_PATH)

    @beartype
    def show(
        self,
        pockets: Optional[list[Pocket]] = None,
        sdf_file: Optional[str] = None,
    ):
        """Visualize the protein structure in a Jupyter notebook using MolStar viewer.

        This method provides interactive 3D visualization of the protein structure with optional
        highlighting of binding pockets and docked ligands. The visualization is rendered directly
        in Jupyter notebooks using the MolStar viewer.

        Args:
            pockets (Optional[list[Pocket]], optional): List of Pocket objects to highlight
                in the visualization. Each pocket will be displayed with its defined color
                and transparency. Defaults to None.
            sdf_file (Optional[str], optional): Path to an SDF file containing docked ligand
                structures. When provided, the ligands will be displayed alongside the protein
                structure. Defaults to None.

        Examples:
            # Visualize protein structure only
            >>> protein = Protein(file="protein.pdb")
            >>> protein.show()

            # Visualize protein with highlighted pockets
            >>> pockets = protein.find_pockets(pocket_count=3)
            >>> protein.show(pockets=pockets)

            # Visualize protein with docked ligands
            >>> protein.show(sdf_file="docked_ligands.sdf")

            # Visualize protein with both pockets and docked ligands
            >>> protein.show(pockets=pockets, sdf_file="docked_ligands.sdf")

        Notes:
            - When pockets are provided, they are displayed with semi-transparent surfaces
              (alpha=0.7) while the protein is shown with a more transparent surface (alpha=0.1)
            - The protein is displayed in cartoon representation when pockets are shown
            - When an SDF file is provided, the visualization includes both the protein and
              the docked ligands in their respective binding poses
        """
        from deeporigin_molstar import JupyterViewer

        current_protein_file = self._dump_state()

        if pockets is None and sdf_file is None:
            protein_viewer = ProteinViewer(
                data=current_protein_file,
                format="pdb",
            )
            html_content = protein_viewer.render_protein()
            JupyterViewer.visualize(html_content)
        elif pockets is not None and sdf_file is None:
            pocket_surface_alpha: float = 0.7
            protein_surface_alpha: float = 0.1

            protein_viewer = ProteinViewer(data=current_protein_file, format="pdb")
            pocket_paths = [str(pocket.file_path) for pocket in pockets]

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
            JupyterViewer.visualize(html_content)
        elif sdf_file is not None:
            from deeporigin_molstar import DockingViewer

            docking_viewer = DockingViewer()
            html_content = docking_viewer.render_with_seperate_crystal(
                protein_data=current_protein_file,
                protein_format="pdb",
                ligands_data=[sdf_file],
                ligand_format="sdf",
            )
            JupyterViewer.visualize(html_content)

    def _repr_html_(self):
        """
        Return the HTML representation of the object for Jupyter Notebook.

        Returns:
            str: The HTML content.
        """
        try:
            if self.info:
                return generate_html_output(self.info)
            return self.visualize()
        except Exception:
            return self.__str__()

    def __str__(self):
        info_str = f"Name: {self.name}\nFile Path: {self.file_path}\n"
        if self.info:
            info_str += f"Info: {self.info}\n"
        return f"Protein:\n  {info_str}"

    def update_coordinates(self, coords: np.ndarray):
        """update coordinates of the protein structure"""

        self.structure.coord = coords

    def get_center_by_residues(self, residues: list[str]) -> np.ndarray:
        """
        Get the center of the protein structure based on specific residues.

        Args:
            residues (list[str]): List of residue names to include in the calculation.
        """
        if not (1 <= len(residues) <= 3):
            print("Please provide 1-3 residue IDs")
            raise ValueError("Invalid number of residue IDs")

        for res_id in residues:
            if not isinstance(res_id, int):
                raise ValueError(f"Residue IDs must be integers. Got: {res_id}")

        mask = np.isin(self.structure.res_id, residues)
        pocket_atoms = self.structure[mask]
        if len(pocket_atoms) == 0:
            raise ValueError(
                f"No atoms found for the specified residue IDs: {residues}"
            )

        warning = ""
        missing_residue_ids = set(residues) - set(pocket_atoms.res_id)
        if missing_residue_ids:
            warning = f"Residue IDs {missing_residue_ids} not found in the structure"

        res_name_id_mapping = {}
        for atom in pocket_atoms:
            res_name_id_mapping[atom.res_name] = atom.res_id

        center = centroid(pocket_atoms)

        with tempfile.TemporaryDirectory() as temp_dir:
            protein_format = "pdb"
            protein_path = os.path.join(temp_dir, "protein.pdb")
            self.to_pdb(protein_path)

            docking_viewer = DockingViewer()
            html = docking_viewer.render_highligh_residues(
                protein_data=protein_path,
                protein_format=protein_format,
                residue_ids=residues,
            )

            if "ATOM" not in html:
                html = ""

        return list(center), warning, JupyterViewer.visualize(html)
