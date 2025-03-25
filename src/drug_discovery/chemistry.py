"""Contains classes and functions for working with molecules, proteins, and related files.

Defines `Ligand` as `Protein` classes, as well as functions for reading/writing SDF files,
SMILES / SDF Conversion, validating data,  `DataFrame` integration, and preparing visualizations.
These can be used together with the `drug_discovery` module for tasks such as docking.

* [`Ligand`][src.drug_discovery.chemistry.Ligand]: Represents a small molecule ligand, accepting either a file path (SDF) or a SMILES
string. Provides `show()` method to display it.
* [`Protein`][src.drug_discovery.chemistry.Protein]: Represents a protein, accepting a local file path (PDB) or a PDB ID.
Provides `show()` method to display it.
"""

import base64
import hashlib
import os
import re
from dataclasses import dataclass, fields
from pathlib import Path
from typing import Optional, Tuple

import pandas as pd
from beartype import beartype
from deeporigin.exceptions import DeepOriginException
from rdkit import Chem


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
            self.properties = read_sdf_properties(self.file)

            # check that there's only one molecule here
            if count_molecules_in_sdf_file(self.file) > 1:
                raise ValueError(
                    "Too many molecules. Expected a single molecule in the SDF file, but got multiple"
                )

        if self.smiles_string is None:
            smiles_string = sdf_to_smiles(self.file)
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

            show_molecules_in_sdf_file(self.file)
        else:
            # only SMILES. use a 2D viewer
            img = smiles_to_base64_png(self.smiles_string)

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


@beartype
def show_molecules_in_sdf_file(sdf_file: str | Path):
    """show molecules in an SDF file in a Jupyter notebook using molstar"""

    from deeporigin_molstar import JupyterViewer, MoleculeViewer

    molecule_viewer = MoleculeViewer(
        data=str(sdf_file),
        format="sdf",
    )
    html_content = molecule_viewer.render_ligand()
    JupyterViewer.visualize(html_content)


@beartype
def show_molecules_in_sdf_files(sdf_files: list[str]):
    """show molecules in an SDF file in a Jupyter notebook using molstar"""

    import tempfile

    temp_dir = tempfile.TemporaryDirectory()

    sdf_file = os.path.join(temp_dir.name, "temp.sdf")

    # combine the SDF files
    merge_sdf_files(sdf_files, sdf_file)

    from deeporigin_molstar import JupyterViewer, MoleculeViewer

    molecule_viewer = MoleculeViewer(
        data=str(sdf_file),
        format="sdf",
    )
    html_content = molecule_viewer.render_ligand()
    JupyterViewer.visualize(html_content)


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
            self.file = download_protein(self.pdb_id)

        self.file = Path(self.file)
        if self.name is None:
            self.name = self.file.name

    def show(self):
        """visualize the protein in a Jupyter notebook using molstar"""

        from deeporigin_molstar import JupyterViewer
        from deeporigin_molstar.src.viewers.protein_viewer import ProteinViewer

        protein_viewer = ProteinViewer(
            data=str(self.file),
            format="pdb",
        )
        html_content = protein_viewer.render_protein()

        JupyterViewer.visualize(html_content)


@beartype
def read_molecules_in_sdf_file(sdf_file: str | Path) -> list[dict]:
    """
    Reads an SDF file containing one or more molecules, and for each molecule:
    - Extracts the SMILES string
    - Extracts all user-defined properties

    Returns:
        list[dict]: A list of dictionaries, where each dictionary has:
            - "smiles_string": str
            - "properties": dict
    """
    from rdkit import Chem

    # Ensure we have a string (RDKit needs a string path)
    sdf_file = str(sdf_file)

    suppl = Chem.SDMolSupplier(sdf_file, sanitize=False)
    if not suppl:  # If the supplier is None or empty
        return []

    output = []
    for mol in suppl:
        if mol is None:
            # Invalid molecule entry in the SDF (could raise a warning or skip)
            continue

        smiles_str = Chem.MolToSmiles(mol)
        properties = {prop: mol.GetProp(prop) for prop in mol.GetPropNames()}

        output.append({"smiles_string": smiles_str, "properties": properties})

    return output


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


def show_ligands(ligands: list[Ligand]):
    """show ligands in the FEP object in a dataframe. This function visualizes the ligands using core-aligned 2D visualizations.

    Args:
        ligands (list[Ligand]): list of ligands

    """

    df = ligands_to_dataframe(ligands)

    # convert SMILES to aligned images
    images = smiles_list_to_base64_png_list(df["Ligand"].tolist())
    df["Ligand"] = images

    # Use escape=False to allow the <img> tags to render as images
    from IPython.display import HTML, display

    display(HTML(df.to_html(escape=False)))


@beartype
def read_sdf_properties(sdf_file: str | Path) -> dict:
    """Reads all user-defined properties from an SDF file (single molecule) and returns them as a dictionary.

    Args:
        sdf_file: Path to the SDF file.

    """

    from rdkit import Chem

    supplier = Chem.SDMolSupplier(
        str(sdf_file),
        sanitize=False,
    )
    mol = supplier[0]  # Assuming a single molecule

    if mol is None:
        raise ValueError("Invalid SDF file or molecule could not be read.")

    return {prop: mol.GetProp(prop) for prop in mol.GetPropNames()}


@beartype
def get_properties_in_sdf_file(sdf_file: str | Path) -> list:
    """Returns a list of all user-defined properties in an SDF file

    Args:
        sdf_file: Path to the SDF file.

    Returns:
        list: A list of the names of all user-defined properties in the SDF file.


    """
    from rdkit import Chem

    # Load molecules from the SDF file
    supplier = Chem.SDMolSupplier(str(sdf_file), sanitize=False)

    properties = []

    for i, mol in enumerate(supplier):
        if mol is None:
            continue  # Skip invalid molecules

        for prop_name in mol.GetPropNames():
            properties.append(prop_name)

    return list(set(properties))


@beartype
def count_molecules_in_sdf_file(sdf_file: str | Path) -> int:
    """
    Count the number of valid (sanitizable) molecules in an SDF file using RDKit,
    while suppressing RDKit's error logging for sanitization issues.

    Args:
        sdf_file: Path to the SDF file.

    Returns:
        int: The number of molecules successfully read in the SDF file.
    """

    from rdkit import Chem, RDLogger

    # Disable RDKit error logging to suppress messages about kekulization/sanitization.
    RDLogger.DisableLog("rdApp.error")

    # Use sanitize=False to defer sanitization until we can handle exceptions ourselves.
    supplier = Chem.SDMolSupplier(str(sdf_file), sanitize=False)
    valid_count = 0
    for mol in supplier:
        if mol is None:
            continue
        try:
            # Manually sanitize the molecule.
            Chem.SanitizeMol(mol)
            valid_count += 1
        except Exception:
            # If sanitization fails, skip this molecule.
            continue
    return valid_count


@beartype
def read_property_values(sdf_file: str | Path, key: str):
    """Given a SDF file with more than 1 molecule, return the values of the properties for each molecule

    Args:
        sdf_file: Path to the SDF file.
        key: The key of the property to read.


    """
    from rdkit import Chem

    suppl = Chem.SDMolSupplier(
        str(sdf_file),
        removeHs=False,
        sanitize=False,
    )
    values = []
    for i, mol in enumerate(suppl, start=1):
        if mol is None:
            value = None
        else:
            if mol.HasProp(key):
                value = mol.GetProp(key).strip()
            else:
                value = None
        values.append(value)
    return values


@beartype
def split_sdf_file(
    *,
    input_sdf_path: str | Path,
    output_prefix: str = "ligand",
    output_dir: Optional[str | Path] = None,
    name_by_property: str = "_Name",
) -> list[Path]:
    """
    Splits a multi-ligand SDF file into individual SDF files, optionally placing
    the output in a user-specified directory. Each output SDF is named using
    the molecule's name (if present) or a fallback prefix.

    Args:
        input_sdf_path: Path to the input SDF file containing multiple ligands.
        output_prefix: Prefix for unnamed ligands. Defaults to "ligand".
        output_dir: Directory to write the output SDF files to. If None,
            output files are written to the same directory as input_sdf_path.

    Returns:
        list[Path]: A list of paths to the generated SDF files.
    """

    from rdkit import Chem

    values = read_property_values(input_sdf_path, name_by_property)
    n_mols = count_molecules_in_sdf_file(input_sdf_path)

    if len(set(values)) != n_mols:
        raise ValueError(
            f"The number of molecules in the SDF file ({n_mols}) is not consistent with the number of values ({len(set(values))}) extracted from the property {name_by_property}. Use a different property to fix this."
        )

    if not isinstance(input_sdf_path, Path):
        input_sdf_path = Path(input_sdf_path)

    if output_dir is None:
        output_dir = input_sdf_path.parent
    else:
        if not isinstance(output_dir, Path):
            output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    suppl = Chem.SDMolSupplier(
        str(input_sdf_path),
        removeHs=False,
        sanitize=False,
    )

    generated_paths = []

    for i, mol in enumerate(suppl, start=1):
        if mol is None:
            continue

        if mol.HasProp(name_by_property):
            mol_name = mol.GetProp(name_by_property).strip()
        else:
            mol_name = f"{output_prefix}_{i}"

        # Replace unsafe filename characters
        safe_name = re.sub(r"[^a-zA-Z0-9_\-]+", "_", mol_name)

        output_file = output_dir / f"{safe_name}.sdf"

        writer = Chem.SDWriter(str(output_file))
        writer.write(mol)
        writer.close()

        generated_paths.append(output_file)

    return generated_paths


@beartype
def smiles_list_to_base64_png_list(
    smiles_list: list[str],
    *,
    size: Tuple[int, int] = (300, 100),
    scale_factor: int = 2,
    reference_smiles: Optional[str] = None,
) -> list[str]:
    """
    Convert a list of SMILES strings to a list of base64-encoded PNG <img> tags.

    This aligns images so that they have consistent core orientation.

    Args:
        smiles_list: List of SMILES strings.
        size: (width, height) of the final rendered image in pixels (CSS downscaled).
        scale_factor: Factor to generate higher-resolution images internally.
        reference_smiles: If provided, all molecules will be oriented to match the 2D layout of this reference molecule.

    """

    from rdkit import Chem
    from rdkit.Chem import AllChem, rdDepictor
    from rdkit.Chem.Draw import rdMolDraw2D

    if reference_smiles is None:
        reference_smiles = smiles_list[0]

    # Prepare the reference molecule
    ref_mol = Chem.MolFromSmiles(reference_smiles)
    if ref_mol is not None:
        AllChem.Compute2DCoords(ref_mol)

    imgs = []

    for smi in smiles_list:
        if not smi:
            imgs.append("N/A")
            continue

        mol = Chem.MolFromSmiles(smi)
        if mol is None:
            imgs.append("N/A")
            continue

        # Generate initial 2D coordinates
        AllChem.Compute2DCoords(mol)

        # If we have a reference molecule, match orientation
        if ref_mol is not None:
            try:
                rdDepictor.GenerateDepictionMatching2DStructure(mol, ref_mol)
            except Exception:
                # In case it fails (incompatible substructures, etc.)
                pass

        # Prepare for drawing
        rdMolDraw2D.PrepareMolForDrawing(mol)

        # Create high-resolution image
        width, height = size[0] * scale_factor, size[1] * scale_factor
        drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
        drawer.DrawMolecule(mol)
        drawer.FinishDrawing()

        png_bytes = drawer.GetDrawingText()
        encoded = base64.b64encode(png_bytes).decode("ascii")

        # Create an <img> tag with CSS-downscaled dimensions
        img_tag = (
            f"<img src='data:image/png;base64,{encoded}' "
            f"style='width:{size[0]}px; height:{size[1]}px;'/>"
        )
        imgs.append(img_tag)

    return imgs


@beartype
def smiles_to_base64_png(
    smiles: str,
    *,
    size=(300, 100),
    scale_factor: int = 2,
) -> str:
    """Convert a SMILES string to an inline base64 <img> tag. Use this if you want to convert a single molecule into an image. If you want to convert a set of SMILES strings (corresponding to a set of related molecules) to images, use `smiles_list_to_base64_png_list`.

    Args:
        smiles (str): SMILES string.
        size (Tuple[int, int], optional): (width, height) of the final rendered image in pixels (CSS downscaled).
        scale_factor (int, optional): Factor to generate higher-resolution images internally.

    """

    from rdkit import Chem
    from rdkit.Chem.Draw import rdMolDraw2D

    if not smiles:
        return "N/A"

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "N/A"

    # Compute high-resolution size
    width, height = size[0] * scale_factor, size[1] * scale_factor

    # Create a high-resolution drawer
    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()

    png = drawer.GetDrawingText()
    encoded = base64.b64encode(png).decode("ascii")

    # Downscale in CSS
    return (
        f"<img src='data:image/png;base64,{encoded}' "
        f"style='width:{size[0]}px; height:{size[1]}px;'/>"
    )


@beartype
def smiles_to_sdf(smiles: str, sdf_path: str) -> None:
    """convert a SMILES string to a SDF file

    Args:
        smiles (str): SMILES string
        sdf_path (str): Path to the SDF file

    """

    from rdkit import Chem
    from rdkit.Chem import AllChem, SDWriter

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Invalid SMILES: {smiles}")

    try:
        Chem.Kekulize(mol)
    except ValueError:
        print(f"Failed to kekulize: {smiles}")

    mol = Chem.AddHs(mol)

    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)

    with SDWriter(sdf_path) as writer:
        writer.write(mol)


@beartype
def sdf_to_smiles(sdf_file: str | Path) -> list[str]:
    """
    Extracts the SMILES strings of all valid molecules from an SDF file using RDKit.

    Args:
        sdf_file (str | Path): Path to the SDF file.

    Returns:
        list[str]: A list of SMILES strings for all valid molecules in the file.
    """
    from rdkit import Chem

    if isinstance(sdf_file, Path):
        sdf_file = str(sdf_file)

    suppl = Chem.SDMolSupplier(sdf_file, sanitize=False)
    if not suppl:
        return []

    smiles_list = []
    for mol in suppl:
        if mol is not None:
            smiles_list.append(Chem.MolToSmiles(mol, canonical=True))

    smiles_list = sorted(set(smiles_list))

    return smiles_list


@beartype
def download_protein(
    pdb_id: str,
    save_dir: str = ".",
) -> str:
    """
    Downloads a PDB structure by its PDB ID from RCSB and saves it to the specified directory.

    Args:
        pdb_id (str): PDB ID of the protein.
        save_dir (str): Directory to save the downloaded PDB file.

    Returns:
        str: Path to the downloaded PDB file.

    Raises:
        Exception: If the download fails.
    """

    from biotite.database.rcsb import fetch

    pdb_id = pdb_id.lower()
    save_dir_path = Path(save_dir)
    save_dir_path.mkdir(parents=True, exist_ok=True)

    file_path = save_dir_path / f"{pdb_id}.pdb"
    if not file_path.exists():
        try:
            fetch(pdb_id, "pdb", save_dir_path)
        except Exception as e:
            raise DeepOriginException(f"Failed to download PDB {pdb_id}: {str(e)}")

    return str(file_path)


@beartype
def merge_sdf_files(
    sdf_file_list: list[str],
    output_path: Optional[str] = None,
) -> str:
    """
    Merge a list of SDF files into a single SDF file.

    Args:
        sdf_file_list (list of str): List of paths to SDF files.

    Returns:
        str: Path to the merged SDF file.
    """
    from rdkit import Chem

    # Get the absolute directory of the first file.
    base_dir = os.path.dirname(os.path.abspath(sdf_file_list[0]))

    # Check that all files are in the same directory.
    for sdf_file in sdf_file_list:
        if os.path.dirname(os.path.abspath(sdf_file)) != base_dir:
            raise ValueError("All input files must be in the same directory.")

    # Create a combined string from the sorted basenames of the input files.
    basenames = sorted([os.path.basename(file) for file in sdf_file_list])
    combined_string = "".join(basenames)

    # Hash the combined string using SHA256 and take the first 10 characters.
    hash_digest = hashlib.sha256(combined_string.encode("utf-8")).hexdigest()[:10]
    output_filename = f"{hash_digest}.sdf"

    if output_path is None:
        output_path = os.path.join(base_dir, output_filename)

    # Check if the output file already exists; if so, do nothing and return it.
    if os.path.exists(output_path):
        return output_path

    # Merge the molecules from all SDF files into the new file.
    writer = Chem.SDWriter(output_path)
    for sdf_file in sdf_file_list:
        supplier = Chem.SDMolSupplier(sdf_file)
        for mol in supplier:
            if mol is not None:  # Skip molecules that failed to parse.
                writer.write(mol)
    writer.close()

    return output_path


@beartype
def filter_sdf_by_smiles(
    *,
    input_sdf_file: str | Path,
    output_sdf_file: str | Path,
    keep_only_smiles: list[str] | pd.Series,
):
    """
    Extracts the SMILES strings of all valid molecules from an SDF file using RDKit.

    Args:
        input_sdf_file (str | Path): Path to the SDF file.
        output_sdf_file (str | Path): Path to the output SDF file.
        keep_only_smiles (list[str] | pd.Series): List or Series of SMILES strings to keep.

    """
    writer = Chem.SDWriter(str(output_sdf_file))

    # Convert pandas Series to list if necessary
    if isinstance(keep_only_smiles, pd.Series):
        keep_only_smiles = keep_only_smiles.tolist()

    keep_only_smiles = [canonicalize_smiles(smiles) for smiles in keep_only_smiles]

    if isinstance(input_sdf_file, Path):
        input_sdf_file = str(input_sdf_file)

    suppl = Chem.SDMolSupplier(
        input_sdf_file,
        sanitize=False,
    )

    for mol in suppl:
        if mol is not None:
            this_smiles = canonicalize_smiles(Chem.MolToSmiles(mol))
            if this_smiles in keep_only_smiles:
                writer.write(mol)

    writer.close()


@beartype
def canonicalize_smiles(smiles: str) -> str:
    """Canonicalize a SMILES string.

    Args:
        smiles (str): SMILES string.

    Returns:
        str: Canonicalized SMILES string.
    """
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"Invalid SMILES: {smiles}")
    return Chem.MolToSmiles(mol, canonical=True)
