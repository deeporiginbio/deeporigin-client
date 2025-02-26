"""module that contains some utility functions for chemistry"""

import base64
import importlib.util
import io
import re
from functools import wraps
from pathlib import Path
from typing import Optional, Union

from beartype import beartype


def requires_rdkit(func):
    """
    A decorator that checks for the presence of RDKit via importlib.util.find_spec.
    If RDKit is unavailable, raises a user-friendly ImportError.
    """

    @wraps(func)
    def wrapper(*args, **kwargs):
        if importlib.util.find_spec("rdkit") is None:
            raise ImportError(
                "RDKit is required for this functionality.\n"
                "Please install it manually \n"
                "or install this package with the extra [tools], for example:\n\n"
                "   pip install deeporigin[tools]\n"
            )
        return func(*args, **kwargs)

    return wrapper


@beartype
@requires_rdkit
def split_sdf_file(
    input_sdf_path: Union[str, Path],
    output_prefix: str = "ligand",
    output_dir: Optional[Union[str, Path]] = None,
) -> None:
    """
    Splits a multi-ligand SDF file into individual SDF files, optionally placing
    the output in a user-specified directory. Each output SDF is named using
    the molecule's name (if present) or a fallback prefix.

    Args:
        input_sdf_path: Path to the input SDF file containing multiple ligands.
        output_prefix: Prefix for unnamed ligands. Defaults to "ligand".
        output_dir: Directory to write the output SDF files to. If None, output files are written to the same directory as input_sdf_path.


    """

    from rdkit import Chem

    # Convert inputs to Path objects if they are not already
    if not isinstance(input_sdf_path, Path):
        input_sdf_path = Path(input_sdf_path)

    if output_dir is None:
        # Default to the directory containing the input file
        output_dir = input_sdf_path.parent
    else:
        if not isinstance(output_dir, Path):
            output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

    # Read the molecules using a string path
    suppl = Chem.SDMolSupplier(str(input_sdf_path), removeHs=False)

    for i, mol in enumerate(suppl, start=1):
        # Skip invalid molecules (i.e., parse failures)
        if mol is None:
            continue

        # Check if the molecule has a name stored in "_Name"
        if mol.HasProp("_Name"):
            mol_name = mol.GetProp("_Name").strip()
        else:
            # Fallback if no name is found
            mol_name = f"{output_prefix}_{i}"

        # Replace unsafe filename characters
        safe_name = re.sub(r"[^a-zA-Z0-9_\-]+", "_", mol_name)

        # Construct the path for the output file
        output_file = output_dir / f"{safe_name}.sdf"

        # RDKit’s writer expects a string path
        writer = Chem.SDWriter(str(output_file))
        writer.write(mol)
        writer.close()

    print(f"Splitting complete! Files written to: {output_dir}")


@beartype
@requires_rdkit
def smiles_to_base64_png(
    smiles: Optional[str],
    size=(300, 10),
) -> str:
    """Convert a SMILES string to an inline base64 <img> tag."""

    from rdkit import Chem
    from rdkit.Chem import Draw

    if not smiles:
        return "N/A"

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        return "N/A"

    img = Draw.MolToImage(mol, size=size, fitImage=False)

    buffer = io.BytesIO()
    img.save(buffer, format="PNG")
    buffer.seek(0)
    encoded = base64.b64encode(buffer.read()).decode("ascii")

    return (
        f"<img src='data:image/png;base64,{encoded}' "
        f"style='max-width:{size[0]}px; max-height:{size[1]}px;'/>"
    )


@beartype
@requires_rdkit
def smiles_to_sdf(smiles: str, sdf_path: str) -> None:
    """convert a SMILES string to a SDF file"""

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
@requires_rdkit
def sdf_to_smiles(sdf_file: str) -> Optional[str]:
    """Extracts the first molecule's SMILES string from an SDF file using RDKit.

    Args:
        sdf_file (str): Path to the SDF file.

    Returns:
        Optional[str]: SMILES string if extraction is successful, else None.
    """

    from rdkit import Chem

    suppl = Chem.SDMolSupplier(sdf_file)
    if not suppl or suppl[0] is None:
        return None

    return Chem.MolToSmiles(suppl[0])


@beartype
@requires_rdkit
def is_ligand_protonated(sdf_file: str) -> bool:
    """
    Determine if the ligand (loaded from an SDF file) is protonated.

    This heuristic function checks for:
      - A positive overall formal charge.
      - The presence of a protonated carboxyl group (-COOH).
      - The presence of a protonated amine (e.g., -NH3+).
      - (Optionally) Other protonated functional groups.

    Parameters:
        sdf_file (str): Path to the SDF file containing the ligand.

    Returns:
        bool: True if protonated evidence is found; False otherwise.
    """
    from rdkit import Chem

    # Load the molecule from the SDF file.
    # Use removeHs=False so that explicit hydrogens are preserved.
    mol = Chem.MolFromMolFile(sdf_file, removeHs=False)
    if mol is None:
        raise ValueError(f"Could not load molecule from {sdf_file}")

    # 1. Check overall formal charge.
    #    (A positive net charge is a strong indicator of protonation.)
    if Chem.GetFormalCharge(mol) > 0:
        return True

    # 2. Check for a protonated carboxylic acid.
    #    The SMARTS '[CX3](=O)[OX1H]' matches a typical protonated -COOH group.
    prot_carboxy_smarts = "[CX3](=O)[OX1H]"
    prot_carboxy = Chem.MolFromSmarts(prot_carboxy_smarts)
    if mol.HasSubstructMatch(prot_carboxy):
        return True

    # 3. Check for a protonated amine.
    #    The SMARTS '[NX3+;H3]' matches, for example, an ammonium group (R-NH3+).
    prot_amine_smarts = "[NX3+;H3]"
    prot_amine = Chem.MolFromSmarts(prot_amine_smarts)
    if mol.HasSubstructMatch(prot_amine):
        return True

    # 4. (Optional) Add more SMARTS queries for other protonated groups
    #    For example, if you want to catch protonated imidazole (imidazolium):
    prot_imidazolium_smarts = "c1[n+][cH]cn1"
    prot_imidazolium = Chem.MolFromSmarts(prot_imidazolium_smarts)
    if mol.HasSubstructMatch(prot_imidazolium):
        return True

    # If no indicators of protonation were found, assume the ligand is not protonated.
    return False
