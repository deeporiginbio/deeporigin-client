"""
LigandModule

Description: This module encapsulates the Ligand class, which represents and manages small molecules (ligands)
in computational biology workflows. It provides functionalities to initialize ligands from various sources,
handle molecular properties, predict ADMET (Absorption, Distribution, Metabolism, Excretion, and Toxicity)
properties, perform protonation state predictions, and visualize ligands within protein structures. The
primary features and methods include:

- **Ligand Class**: Represents a ligand molecule, managing its data loaded from identifiers, file paths,
  SMILES strings, or direct block content. Handles molecule creation, property assignment, and ensures data
  integrity during initialization.

- **Initialization**: Allows creating a Ligand instance using an identifier (e.g., PubChem ID), local file,
  SMILES string, or block content. Ensures only one source is provided and properly parses the molecular data.

- **Property Management**: Stores and manages ligand-specific properties such as volume, SASA (Solvent Accessible
  Surface Area), hydrophobicity, polarity, and drugability scores. Provides methods to set and retrieve these
  properties, facilitating comprehensive molecular analysis.

- **ADMET Prediction**: Integrates with predictive models to estimate ADMET properties of ligands, aiding in the
  assessment of their drug-likeness and suitability for therapeutic applications.

- **Protonation State Prediction**: Predicts the protonation states of ligands at specified pH levels, providing
  insights into their chemical behavior and interactions within biological environments.

- **Visualization**: Integrates with the MoleculeViewer to render interactive visualizations of ligands within
  Jupyter Notebooks, facilitating intuitive analysis and presentation of molecular structures.

- **Utility Methods**: Includes static and class methods for creating Ligand instances from SDF and CSV files,
  converting molecular formats, fetching SMILES strings from APIs, and handling bulk visualization of multiple
  ligands.

Dependencies: Utilizes libraries such as RDKit for cheminformatics, Pandas for data manipulation, Termcolor for
colored terminal output, Requests for API interactions, and integrates with the MoleculeViewer from the
deeporigin_molstar package for visualization purposes.

Usage Example:
```python
# Initialize a Ligand instance from a SMILES string
ligand = Ligand(smiles="CCO", name="Ethanol")

# Predict ADMET properties for the ligand
admet = ligand.admet_properties()

# Predict protonation states at pH 7.4
ligand.protonate(pH=7.4)

# Display ligand properties
print(ligand)

# Visualize the ligand within a protein structure
ligand.visualize()
"""

from dataclasses import dataclass, field
import os
from pathlib import Path
import tempfile
from typing import Any, List, Optional, Union

from beartype import beartype
from deeporigin_molstar import MoleculeViewer
import numpy as np
import pandas as pd
from rdkit import Chem
from tqdm import tqdm

from deeporigin.drug_discovery.structures.internal_structures import (
    Molecule,
    mol_from_block,
    mol_from_smiles,
)
from deeporigin.drug_discovery.utilities.visualize import jupyter_visualization
from deeporigin.exceptions import DeepOriginException

from .entity import Entity


@dataclass
class Ligand(Entity):
    """A class representing a ligand molecule in drug discovery workflows.
    The Ligand class provides functionality to create, manipulate, and analyze small molecules
    (ligands) in computational drug discovery. It supports various input formats and provides
    methods for property prediction, visualization, and file operations.

    Attributes:
        identifier (Optional[str]): Ligand identifier (e.g., PubChem ID)
        file_path (Optional[str]): Path to the ligand file
        smiles (Optional[str]): SMILES string representing the ligand
        block_type (Optional[str]): Format of the block content ('mol', 'mol2', 'sdf', 'pdb')
        block_content (Optional[str]): String containing the molecule data
        name (Optional[str]): Optional name of the ligand
        seed (Optional[int]): Random seed for coordinate generation
        xref_protein (Optional[str]): Cross-reference to protein
        xref_ins_code (Optional[str]): Cross-reference insertion code
        xref_residue_id (Optional[str]): Cross-reference residue ID
        xref_protein_chain_id (Optional[str]): Cross-reference protein chain ID
        save_to_file (bool): Whether to save the ligand to file
        properties (dict): Dictionary of ligand properties
        mol (Optional[Molecule]): Direct Molecule object initialization


    Examples:
        >>> # Create from SMILES
        >>> ligand = Ligand.from_smiles("CCO", name="Ethanol")

        >>> # Create from SDF file
        >>> ligand = Ligand.from_sdf("ligand.sdf")

        >>> # Get properties
        >>> center = ligand.get_center()
        >>> props = ligand.admet_properties()

        >>> # Visualize
        >>> ligand.visualize()

        >>> # Save to file
        >>> ligand.write_to_file("output.pdb")
    """

    identifier: str | None = None
    file_path: str | None = None
    smiles: str | None = None
    block_type: str | None = None
    block_content: str | None = None
    name: str | None = None
    seed: int | None = None
    xref_protein: str | None = None
    xref_ins_code: str | None = None
    xref_residue_id: str | None = None
    xref_protein_chain_id: str | None = None
    save_to_file: bool = False
    properties: dict = field(default_factory=dict)
    mol: Molecule | None = None

    # Additional attributes that are initialized in __post_init__
    protonated_smiles: str | None = field(init=False, default=None)
    hac: int = field(init=False, default=0)
    available_for_docking: bool = field(init=False, default=True)

    _remote_path_base = "entities/ligands/"

    @classmethod
    @beartype
    def from_rdkit_mol(
        cls,
        mol: Chem.rdchem.Mol,
        name: str = "",
        save_to_file: bool = False,
        **kwargs: Any,
    ) -> "Ligand":
        """
        Create a Ligand instance from an RDKit Mol object.

        Args:
            mol (Chem.rdchem.Mol): RDKit molecule object to convert to a Ligand
            name (str, optional): Name of the ligand. Defaults to "".
            save_to_file (bool, optional): Whether to save the ligand to file. Defaults to False.
            **kwargs: Additional arguments to pass to the constructor

        Returns:
            Ligand: A new Ligand instance initialized from the RDKit molecule

        Example:
            >>> from rdkit import Chem
            >>> mol = Chem.MolFromSmiles("CCO")
            >>> ligand = Ligand.from_rdkit_mol(mol, name="Ethanol")
        """
        # Get name from properties if available
        if mol.HasProp("_Name"):
            name = mol.GetProp("_Name")
        elif name == "" and "properties" in kwargs and "_Name" in kwargs["properties"]:
            name = kwargs["properties"]["_Name"]

        return cls(
            mol=Molecule(mol, name=name),
            name=name,
            save_to_file=save_to_file,
            **kwargs,
        )

    @classmethod
    def from_smiles(
        cls,
        smiles: str,
        name: str = "",
        save_to_file: bool = False,
        **kwargs: Any,
    ) -> "Ligand":
        """
        Create a Ligand instance from a SMILES string.

        Args:
            smiles (str): SMILES string representing the ligand
            name (str, optional): Name of the ligand. Defaults to "".
            save_to_file (bool, optional): Whether to save the ligand to file. Defaults to False.
            **kwargs: Additional arguments to pass to the constructor

        Returns:
            Ligand: A new Ligand instance

        Example:
            >>> ligand = Ligand.from_smiles(
            ...     smiles="CCO",  # Ethanol
            ...     name="Ethanol",
            ...     save_to_file=False
            ... )
            >>> print(ligand.smiles)
            CCO
        """
        try:
            # Create a Molecule object from the SMILES string
            mol = mol_from_smiles(smiles)
        except ValueError as e:
            raise DeepOriginException(str(e)) from e

        return cls(
            mol=mol,
            smiles=smiles,
            name=name,
            save_to_file=save_to_file,
            **kwargs,
        )

    @classmethod
    def from_block_content(
        cls,
        block_content: str,
        block_type: str,
        name: str = "",
        save_to_file: bool = False,
        **kwargs: Any,
    ) -> "Ligand":
        """
        Create a Ligand instance from block content.

        Args:
            block_content (str): String containing the molecule data
            block_type (str): Format of the block content ('mol', 'mol2', 'sdf', 'pdb')
            name (str, optional): Name of the ligand. Defaults to "".
            save_to_file (bool, optional): Whether to save the ligand to file. Defaults to False.
            **kwargs: Additional arguments to pass to the constructor

        Returns:
            Ligand: A new Ligand instance
        """

        mol = mol_from_block(block_type, block_content)

        return cls(
            block_content=block_content,
            block_type=block_type,
            mol=mol,
            name=name,
            save_to_file=save_to_file,
            **kwargs,
        )

    @classmethod
    def from_identifier(
        cls,
        identifier: str,
        name: Optional[str] = None,
        save_to_file: bool = False,
        **kwargs: Any,
    ) -> "Ligand":
        """
        Create a Ligand instance from a chemical identifier.

        Args:
            identifier (str): Chemical identifier (e.g., common name, PubChem name, drug name)
            name (str, optional): Name of the ligand. If not provided, uses the identifier. Defaults to "".
            save_to_file (bool, optional): Whether to save the ligand to file. Defaults to False.
            **kwargs: Additional arguments to pass to the constructor

        Returns:
            Ligand: A new Ligand instance initialized from the chemical identifier

        Example:
            >>> # Create ATP molecule
            >>> atp = Ligand.from_identifier("ATP", name="ATP")
            >>>
            >>> # Create serotonin molecule
            >>> serotonin = Ligand.from_identifier(
            ...     identifier="serotonin",
            ...     name="Serotonin"
            ... )

        Raises:
            DeepOriginException: If the identifier cannot be resolved to a valid molecule
        """
        try:
            mol = Molecule.from_smiles_or_name(
                name=identifier,
                add_coords=True,
            )
        except Exception as e:
            raise DeepOriginException(
                f"Could not resolve chemical identifier '{identifier}': {str(e)}"
            ) from e

        if name is None:
            name = identifier

        return cls(
            mol=mol,
            identifier=identifier,
            name=name,
            save_to_file=save_to_file,
            **kwargs,
        )

    @classmethod
    def from_sdf(
        cls,
        file_path: str,
        *,
        sanitize: bool = True,
        removeHs: bool = False,
    ) -> Union[List["Ligand"], "Ligand"]:
        """
        Create Ligand instances from an SDF file.

        Args:
            file_path (str): The path to the SDF file.

        Returns:
            list[Ligand]: A list of Ligand instances created from the SDF file.

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If the file cannot be parsed correctly.
        """
        path = Path(file_path)
        if not path.exists():
            raise FileNotFoundError(f"The file '{file_path}' does not exist.")

        ligands = []
        try:
            suppl = Chem.SDMolSupplier(str(path), sanitize=sanitize, removeHs=removeHs)
            for idx, mol in enumerate(suppl, start=1):
                try:
                    if mol is None:
                        print(
                            f"Warning: Skipping molecule at index {idx} due to parsing error."
                        )
                        continue
                    ligand = Ligand.from_rdkit_mol(
                        mol,
                        properties=mol.GetPropsAsDict(),
                    )
                    ligands.append(ligand)
                except Exception as e:
                    print(
                        f"Error: Failed to create Ligand from SDF file molecule_idx = '{idx}': {str(e)}"
                    )
        except Exception as e:
            raise ValueError(
                f"Failed to create Ligands from SDF file '{file_path}': {str(e)}"
            ) from e

        if len(ligands) == 1:
            # Set file_path for single ligand case
            ligands[0].file_path = str(path)
            return ligands[0]
        return ligands

    @classmethod
    def from_csv(
        cls,
        file_path: str,
        smiles_column: str = "smiles",
    ) -> list["Ligand"]:
        """
        Create Ligand instances from a CSV file containing SMILES strings and additional properties.

        Args:
            file_path (str): The path to the CSV file.
            smiles_column (str, optional): The name of the column containing SMILES strings. Defaults to "smiles".

        Returns:
            list[Ligand]: A list of Ligand instances created from the CSV file.

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If the CSV does not contain the specified smiles column or if SMILES strings are invalid.
        """
        path = Path(file_path)
        if not path.exists():
            raise FileNotFoundError(f"The file '{file_path}' does not exist.")

        # First read just the header to check for the smiles column
        df_header = pd.read_csv(file_path, nrows=0)
        if smiles_column not in df_header.columns:
            raise ValueError(
                f"Column '{smiles_column}' not found in CSV file '{file_path}'. Available columns: {', '.join(df_header.columns)}"
            )

        ligands = []
        try:
            df = pd.read_csv(file_path)
            normalized_columns = [col.strip().lower() for col in df.columns]

            if smiles_column.lower() not in normalized_columns:
                raise ValueError(f"CSV file must contain a '{smiles_column}' column.")

            smiles_col_index = normalized_columns.index(smiles_column.lower())
            smiles_col = df.columns[smiles_col_index]
            other_columns = [col for col in df.columns if col != smiles_col]

            for idx, row in df.iterrows():
                try:
                    smiles = row[smiles_col]
                    if pd.isna(smiles):
                        print(
                            f"Warning: Skipping row {idx + 1}: SMILES value is missing."
                        )
                        continue
                    mol = Chem.MolFromSmiles(smiles)
                    if mol is None:
                        print(
                            f"Warning: Skipping row {idx + 1}: Invalid SMILES '{smiles}'."
                        )
                        continue

                    # Create properties dictionary
                    properties = {}
                    for col in other_columns:
                        value = row[col]
                        if pd.notna(value):
                            properties[col] = value

                    # Get name from properties if available
                    name = properties.get("Name", "")

                    # Create ligand using from_smiles
                    ligand = cls.from_smiles(
                        smiles=smiles,
                        name=name,
                        properties=properties,
                    )
                    ligands.append(ligand)
                except Exception as e:
                    print(
                        f"Error: Failed to create Ligand from CSV file row {idx + 1}: {str(e)}"
                    )

        except pd.errors.EmptyDataError as e:
            raise ValueError(f"The CSV file '{file_path}' is empty.") from e
        except pd.errors.ParserError as e:
            raise ValueError(f"Error parsing CSV file '{file_path}': {str(e)}") from e
        except Exception as e:
            raise ValueError(
                f"Failed to create Ligands from CSV file '{file_path}': {str(e)}"
            ) from e

        return ligands

    def __post_init__(self):
        """
        Initialize a Ligand instance from an identifier, file path, SMILES string,
        block content, or direct Molecule object.
        """

        # check that a mol exists
        if self.mol is None:
            raise ValueError(
                "mol must be provided when initializing from an identifier, file path, SMILES string, or block content."
            )

        if self.smiles is None:
            self.smiles = self.mol.smiles

        if self.smiles is not None:
            self.mol = mol_from_smiles(self.smiles)
            self.block_type = "mol"
            self.block_content = self.mol.molblock()
        else:
            raise ValueError("No valid source provided for ligand initialization.")

        if self.mol is None:
            raise ValueError("Failed to create molecule.")

        self.smiles = self.mol.smiles

        self.name = self.mol.name if self.mol.name else self.name or "Unknown_Ligand"
        directory = Path(self._get_directory())
        if self.name == "Unknown_Ligand":
            num = len(list(directory.glob(f"{self.name}*")))
            self.name = f"{self.name}_{num + 1}"

        self.hac = self.mol.m.GetNumHeavyAtoms()
        if self.hac < 5:
            print("Warning: Ligand has less than 5 heavy atoms.")

        file_props = self.mol.m.GetPropsAsDict()
        for key, value in file_props.items():
            self.properties[key] = value

        self.available_for_docking = not self.mol.contains_boron
        if self.save_to_file:
            self.write_to_file(output_format="sdf")

    @property
    def coordinates(self):
        return np.array(self.mol.coords(), dtype=np.float32)

    @property
    def atom_types(self):
        return self.mol.species()

    def set_property(self, prop_name: str, prop_value):
        """
        Set a property for the ligand molecule.

        Parameters:
        - prop_name (str): Name of the property.
        - prop_value: Value of the property.

        Example:
        ```python
        ligand.set_property("BindingAffinity", 5.6)
        ```
        """
        self.properties[prop_name] = prop_value
        self.mol.m.SetProp(prop_name, str(prop_value))

    def get_property(self, prop_name: str):
        """
        Get the value of a property for the ligand molecule.

        Parameters:
        - prop_name (str): Name of the property to retrieve.

        Returns:
        - The value of the property if it exists, otherwise None.

        Example:
        ```python
        binding_affinity = ligand.get_property("BindingAffinity")
        ```
        """
        value = self.properties.get(prop_name)
        if value is not None:
            return value

        if self.mol.m.HasProp(prop_name):
            value = self.mol.m.GetProp(prop_name)
            self.properties[prop_name] = value
            return value

        return None

    def write_to_file(self, output_path: str = "", output_format: str = ""):
        """
        Writes the ligand molecule to a file, including all properties.

        Parameters:
        - output_path (str): Path where the ligand will be written.

        Raises:
        - ValueError: If the file extension is unsupported.
        - Exception: If writing to the file fails.

        Example:
        ```python
        ligand.write_to_file('/path/to/output.pdb')
        ```
        """
        try:
            if output_format == "" and output_path == "":
                raise ValueError("Please provide either output_path or output_format.")

            if not output_path:
                output_path = str(
                    Path(self._get_directory()) / f"{self.name}.{output_format}"
                )

            path = Path(output_path)
            extension = path.suffix.lower()
            if not output_format:
                output_format = extension

            if output_format and output_format[0] != ".":
                output_format = f".{output_format}"

            if extension and extension != output_format:
                print(
                    "Warning: Output format does not match the file extension. Writing to provided output format."
                )
                extension = output_format

            if self.name is not None:
                self.set_property("_Name", self.name)
            if self.mol.smiles is not None:
                self.set_property("_SMILES", self.mol.smiles)
            if self.properties:
                for prop_name, prop_value in self.properties.items():
                    self.set_property(prop_name, str(prop_value))

            if extension == ".pdb":
                pdb_block = Chem.MolToPDBBlock(self.mol.m)
                remark_lines = ""
                for prop_name, prop_value in self.mol.m.GetPropsAsDict().items():
                    remark_lines += f"REMARK   {prop_name}: {prop_value}\n"
                pdb_block_with_remarks = remark_lines + pdb_block
                path.write_text(pdb_block_with_remarks)
            elif extension == ".sdf":
                with tempfile.NamedTemporaryFile(
                    mode="w+", suffix=".sdf", delete=False
                ) as temp_file:
                    writer = Chem.SDWriter(temp_file.name)
                    writer.write(self.mol.m)
                    writer.close()
                    temp_file.flush()
                    temp_file.seek(0)
                    path.write_text(temp_file.read())
            elif extension == ".mol":
                mol_block = Chem.MolToMolBlock(self.mol.m)
                prop_lines = ""
                for prop_name, prop_value in self.mol.m.GetPropsAsDict().items():
                    prop_lines += f">  <{prop_name}>\n{prop_value}\n\n"
                mol_block_with_props = mol_block + "\n" + prop_lines
                path.write_text(mol_block_with_props)
            else:
                raise ValueError(
                    f"Unsupported file extension '{extension}'. Supported extensions are '.pdb', '.mol', '.sdf'."
                )

        except Exception as e:
            raise ValueError(
                f"Failed to write structure to file {output_path}: {str(e)}"
            ) from e

    def get_center(self) -> Optional[list[float]]:
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
            print("Warning: Coordinates are not available for this ligand.")
            return None
        center = self.coordinates.mean(axis=0)
        return [float(x) for x in center.tolist()]

    def draw(self):
        """
        Draw the ligand molecule.

        Example:
        ```python
        ligand.draw()
        ```
        """
        return self.mol.draw()

    @jupyter_visualization
    def show(self) -> str:
        """
        Visualize the current state of the ligand molecule.

        Returns:
        - str: HTML representation of the visualization.

        Raises:
        - Exception: If visualization fails.

        Example:
        ```python
        ligand.show()
        ```
        """
        try:
            temp_file = Path(tempfile.gettempdir()) / f"{self.name}_visualize.sdf"
            self.write_to_file(str(temp_file))

            viewer = MoleculeViewer(str(temp_file), format="sdf")
            ligand_config = viewer.get_ligand_visualization_config()
            html = viewer.render_ligand(ligand_config=ligand_config)

            return html
        except Exception as e:
            raise ValueError(f"Visualization failed: {str(e)}") from e

    @classmethod
    @beartype
    def convert_to_sdf(
        cls,
        block_content: str,
        block_type: str,
    ) -> str:
        """
        Convert a ligand block content to SDF format.

        Args:
            block_content (str): The block content of the ligand.
            block_type (str): The type of the block content.

        Returns:
            str: The ligand block content in SDF format.
        """
        try:
            molecule = mol_from_block(
                block_type,
                block_content,
                sanitize=True,
                remove_hs=False,
            )
            with tempfile.TemporaryFile(mode="w+", suffix=".sdf") as temp_file:
                writer = Chem.SDWriter(temp_file.name)
                writer.write(molecule.m)
                writer.close()

            return molecule.molblock()
        except Exception as e:
            raise ValueError(
                f"Failed to convert ligand block content to SDF: {str(e)}"
            ) from e

    @classmethod
    @jupyter_visualization
    def visualize_ligands_from_sdf(cls, file_path: str):
        """
        Visualize ligands from an SDF file.

        Args:
            file_path (str): The path to the SDF file.

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If the file cannot be parsed correctly.
        """
        try:
            viewer = MoleculeViewer(str(file_path), format="sdf")
            ligand_config = viewer.get_ligand_visualization_config()
            html = viewer.render_ligand(ligand_config=ligand_config)

            return html
        except Exception as e:
            raise ValueError(f"Visualization failed: {str(e)}") from e

    @classmethod
    @jupyter_visualization
    def visualize_ligands(cls, ligands: list["Ligand"]):
        """
        Visualize ligands.

        Args:
            ligands: list["Ligand"]: The list of ligands objects to visualize.

        Raises:
            FileNotFoundError: If the file does not exist.
            ValueError: If the file cannot be parsed correctly.
        """
        try:
            sdf_data = []
            current_file = f"{tempfile.mkstemp()[1]}.sdf"
            for ligand in ligands:
                ligand.write_to_file(output_format="sdf", output_path=current_file)

                with open(current_file, "r") as fd:
                    data = fd.read()

                sdf_data.append(data)

            sdf_data = "".join(sdf_data)
            # Write the consolidated SDF data to a temporary file
            with open(current_file, "w") as fd:
                fd.write(sdf_data)

            # Use visualize_ligands_from_sdf to handle the visualization
            return cls.visualize_ligands_from_sdf(current_file)
        except Exception as e:
            raise ValueError(f"Visualization failed: {str(e)}") from e

    def _repr_html_(self) -> str:
        """
        Return the HTML representation of the object for Jupyter Notebook.

        Returns:
            str: The HTML content.
        """
        try:
            print(self.mol.m)
            return self.show()
        except Exception as e:
            print(f"Warning: Failed to generate HTML representation: {str(e)}")
            return self.__str__()

    def __str__(self) -> str:
        info_str = (
            f"Name: {self.name}\nSMILES: {self.mol.smiles}\nHeavy Atoms: {self.hac}\n"
        )
        if self.properties:
            info_str += "Properties:\n"
            for prop_name, prop_value in self.properties.items():
                info_str += f"  {prop_name}: {prop_value}\n"

        if self.xref_protein is not None:
            info_str += (
                f"Cross-reference Protein Chain ID: {self.xref_protein_chain_id}\n"
            )
            info_str += f"Cross-reference Residue ID: {self.xref_residue_id}\n"
            info_str += f"Cross-reference Insertion Code: {self.xref_ins_code}\n"

        return f"Ligand:\n  {info_str}"

    def __repr__(self) -> str:
        return self.__str__()

    @staticmethod
    def _get_directory() -> str:
        """
        Generates and ensures the existence of a directory for ligands.

        Returns:
            str: The path to the ligands directory (~/.deeporigin/ligands).
        """
        ligands_base_dir = Path.home() / ".deeporigin" / "ligands"
        ligands_base_dir.mkdir(parents=True, exist_ok=True)

        return str(ligands_base_dir)

    def admet_properties(self) -> str:
        """
        Predict ADMET properties for the ligand.

        Returns:
            str: A string containing the predicted ADMET properties.
        """

        from deeporigin.functions.molprops import molprops

        try:
            props = molprops(self.mol.smiles)["properties"]
            for key, value in props.items():
                self.set_property(key, value)

            return props
        except Exception as e:
            raise ValueError(f"Failed to predict ADMET properties: {str(e)}") from e

    def protonate(self, pH: float = 7.4, filter_percentage: float = 1):
        """
        Predicts the right protonation of a molecule at given pH value.

        Parameters:
        - entry: A single or multiple ligands represented as SMILES or Ligand instances.
        - pH: pH value of the solvent for concentration calculation. Default is 7.4.
        - filter_percentage: Percentage threshold for filtering low concentration states. Default is 1.

        Returns:
        - ProtonationReport: A ProtonationReport instance.
        """

        raise NotImplementedError("Protonation prediction not implemented yet.")
        # try:
        #     smiles = protonate(
        #         pH=pH,
        #         smiles=self.mol.smiles,
        #         filter_percentage=filter_percentage,
        #     )
        #     if smiles:
        #         self.protonated_smiles = smiles
        #         self.set_property("ProtonatedSMILES", smiles)
        # except Exception as e:
        #     raise ValueError(f"Failed to protonate the ligand molecule: {str(e)}")

    def update_coordinates(self, coords: np.ndarray):
        """update coordinates of the ligand structure"""

        if self.mol.m.GetNumConformers() == 0:
            raise ValueError("Ligand molecule has no conformers to update.")

        conformer = self.mol.m.GetConformer()
        mol_without_hs = Chem.RemoveHs(self.mol.m)

        conformer_no_hs = mol_without_hs.GetConformer()
        if coords.shape[0] != conformer.GetNumAtoms():
            if coords.shape[0] != conformer_no_hs.GetNumAtoms():
                raise ValueError(
                    "Number of ligand atoms does not match the conformer's atom count."
                )

        conformer.SetPositions(coords.astype(np.float64))

    @classmethod
    def protonate_molecules(cls, ligands):
        """
        Predicts the right protonation of a molecule at given pH value.

        Parameters:
        - entry: A single or multiple ligands represented as SMILES or Ligand instances.
        - pH: pH value of the solvent for concentration calculation. Default is 7.4.
        - filter_percentage: Percentage threshold for filtering low concentration states. Default is 1.

        Returns:
        - ProtonationReport: A ProtonationReport instance.
        """
        mols = []

        for i in tqdm(range(0, len(ligands)), desc="Protonating Molecules"):
            ligand = ligands[i]
            if isinstance(ligand, str):
                try:
                    ligand = Ligand(smiles=ligand)
                except Exception as e:
                    print(f"Error: Failed to create Ligand from SMILES: {str(e)}")
                    continue
            try:
                if not ligand.protonated_smiles:
                    ligand.protonate()
            except Exception as e:
                print(f"Error: Failed to protonate the ligand molecule: {str(e)}")
                continue

            mols.append(ligand)
        return mols


@beartype
def ligands_to_dataframe(ligands: list[Ligand]):
    """convert a list of ligands to a pandas dataframe"""

    import pandas as pd

    smiles_list = [ligand.smiles for ligand in ligands]
    file_list = [
        os.path.basename(ligand.file_path) if ligand.file_path is not None else None
        for ligand in ligands
    ]

    data = {
        "Ligand": smiles_list,
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


def show_ligands(ligands: list[Ligand]):
    """show ligands in the FEP object in a dataframe. This function visualizes the ligands using core-aligned 2D visualizations.

    Args:
        ligands (list[Ligand]): list of ligands

    """

    df = ligands_to_dataframe(ligands)

    from rdkit.Chem import PandasTools

    PandasTools.AddMoleculeColumnToFrame(df, smilesCol="Ligand", molCol="Ligand")
    PandasTools.RenderImagesInAllDataFrames()

    # show structure first
    new_order = ["Ligand"] + [col for col in df.columns if col != "Ligand"]

    # re‑index your DataFrame
    df = df[new_order]

    from deeporigin.utils.notebook import get_notebook_environment

    if get_notebook_environment() == "marimo":
        import marimo as mo

        return mo.plain(df)

    else:
        return df
