"""module that contains some utility functions for chemistry"""

import base64
import importlib.util
import re
from dataclasses import dataclass, fields
from functools import wraps
from pathlib import Path
from typing import List, Optional, Tuple, Union

from beartype import beartype


@dataclass
class Ligand:
    """class to represent a ligand (typically backed by a SDF file)"""

    file: Union[str, Path]
    smiles_string: Optional[str] = None
    n_molecules: Optional[int] = None

    # this ID keeps track of whether it is uploaded to deep origin or not
    _do_id: Optional[str] = None

    # this stores user-defined properties
    properties: Optional[dict] = None

    def __post_init__(self):
        """generates a SMILES if it doesn't exist"""

        # read user-defined properties
        self.properties = read_sdf_properties(self.file)

        # check that there's only one molecule here
        if self.n_molecules is None:
            if count_molecules_in_sdf_file(self.file) > 1:
                raise ValueError(
                    "Too many molecules. Expected a single molecule in the SDF file, but got multiple"
                )
            self.n_molecules = 1

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


@dataclass
class Protein:
    """class to represent a protein (typically backed by a PDB file)"""

    file: Union[str, Path]
    name: Optional[str] = None

    # this ID keeps track of whether it is uploaded to deep origin or not
    _do_id: Optional[str] = None

    def __post_init__(self):
        self.file = Path(self.file)
        if self.name is None:
            self.name = self.file.name


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
def read_sdf_properties(sdf_file: str | Path) -> dict:
    """Reads all user-defined properties from an SDF file (single molecule) and returns them as a dictionary."""

    from rdkit import Chem

    supplier = Chem.SDMolSupplier(str(sdf_file), sanitize=False)
    mol = supplier[0]  # Assuming a single molecule

    if mol is None:
        raise ValueError("Invalid SDF file or molecule could not be read.")

    return {prop: mol.GetProp(prop) for prop in mol.GetPropNames()}


@beartype
@requires_rdkit
def get_properties_in_sdf_file(sdf_file: str) -> list:
    """returns a list of all user-defined properties in an SDF file"""
    from rdkit import Chem

    # Load molecules from the SDF file
    supplier = Chem.SDMolSupplier(sdf_file, sanitize=False)

    properties = []

    for i, mol in enumerate(supplier):
        if mol is None:
            continue  # Skip invalid molecules

        for prop_name in mol.GetPropNames():
            properties.append(prop_name)

    return list(set(properties))


@beartype
@requires_rdkit
def count_molecules_in_sdf_file(sdf_file: str | Path) -> int:
    """
    Count the number of valid (sanitizable) molecules in an SDF file using RDKit,
    while suppressing RDKit's error logging for sanitization issues.

    Args:
        sdf_file (str or Path): Path to the SDF file.

    Returns:
        int: The number of molecules successfully sanitized.
    """
    # Disable RDKit error logging to suppress messages about kekulization/sanitization.

    from rdkit import Chem, RDLogger

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
@requires_rdkit
def read_property_values(sdf_file: str, key: str):
    """given a SDF file with more than 1 molecule, return the values of the properties for each molecule"""
    from rdkit import Chem

    suppl = Chem.SDMolSupplier(
        sdf_file,
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
@requires_rdkit
def split_sdf_file(
    input_sdf_path: Union[str, Path],
    output_prefix: str = "ligand",
    output_dir: Optional[Union[str, Path]] = None,
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
        List[Path]: A list of paths to the generated SDF files.
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
@requires_rdkit
def smiles_list_to_base64_png_list(
    smiles_list: List[str],
    *,
    size: Tuple[int, int] = (300, 100),
    scale_factor: int = 2,
    reference_smiles: Optional[str] = None,
) -> List[str]:
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
    ref_mol = None
    if reference_smiles:
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
@requires_rdkit
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
def sdf_to_smiles(sdf_file: Union[str, Path]) -> list[str]:
    """
    Extracts the SMILES strings of all valid molecules from an SDF file using RDKit.

    Args:
        sdf_file (Union[str, Path]): Path to the SDF file.

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
            smiles_list.append(Chem.MolToSmiles(mol))

    return smiles_list


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
