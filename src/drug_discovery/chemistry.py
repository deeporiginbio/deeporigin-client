"""Contains functions for working with SDF files."""

import hashlib
import os
from pathlib import Path
import re
from typing import Optional

from beartype import beartype
from rdkit import Chem
from rdkit.Chem import AllChem


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
    for _, mol in enumerate(suppl, start=1):
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
def smiles_to_sdf(smiles: str, sdf_path: str) -> None:
    """convert a SMILES string to a SDF file

    Args:
        smiles (str): SMILES string
        sdf_path (str): Path to the SDF file

    """

    from rdkit import Chem
    from rdkit.Chem import SDWriter

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
