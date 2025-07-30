"""Contains functions for working with SDF files."""

import base64
import hashlib
import os
from pathlib import Path
import re
from typing import Optional


#from beartype import beartype
from typing import List, Tuple

from rdkit import Chem
from rdkit.Chem import rdFMCS, AllChem
from rdkit.Chem import rdForceFieldHelpers
from rdkit.Chem.Scaffolds import MurckoScaffold
import pandas as pd

def read_molecules(input_file: str, column: str = None, threshold: float = 0.9) -> List[Chem.Mol]:
    """
    Read molecules from CSV/TSV/SDF/MOL/SMILES files into RDKit Mol objects.

    Args:
        input_file (str): Path to a .csv, .tsv, .sdf, .mol, .smi or .smiles file.
        column (Optional[str]): If CSV/TSV, the name of the SMILES column. Auto‑detected if None.
        threshold (float): Fraction of valid SMILES needed to auto‑detect a column.

    Returns:
        List[Chem.Mol]: Parsed RDKit Mol objects.
    """
    ext = os.path.splitext(input_file)[-1].lower()
    mols: List[Chem.Mol] = []

    if ext in [".csv", ".tsv"]:
        sep = "\t" if ext == ".tsv" else ","
        df = pd.read_csv(input_file, sep=sep)

        if column:
            smiles_list = df[column].dropna().astype(str)
            mols = [Chem.MolFromSmiles(smi) for smi in smiles_list]
        else:
            found_column = None
            for col in df.columns:
                smiles_list = df[col].dropna().astype(str)
                parsed = [Chem.MolFromSmiles(smi) for smi in smiles_list]
                valid_count = sum(1 for m in parsed if m is not None)
                if len(smiles_list) > 0 and (valid_count / len(smiles_list) >= threshold):
                    found_column = col
                    mols = [m for m in parsed if m is not None]
                    break
            if found_column is None:
                raise ValueError("Could not find a valid SMILES column in the input file.")
            print(f"Detected SMILES column: '{found_column}'")

    elif ext in [".sdf", ".mol"]:
        suppl = Chem.SDMolSupplier(input_file, removeHs=False)
        mols = [m for m in suppl if m is not None]

    elif ext in [".smi", ".smiles"]:
        with open(input_file, "r") as f:
            lines = f.read().splitlines()
        for line in lines:
            smi = line.strip().split()[0]
            m = Chem.MolFromSmiles(smi)
            if m is not None:
                mols.append(m)

    else:
        raise ValueError("Unsupported file format: use .csv, .tsv, .sdf, .mol, or .smi/.smiles.")

    print(f"Number of parsed molecules: {len(mols)}")
    return mols



#@beartype
def read_molecules_in_sdf_file(sdf_file: str | Path) -> list[dict]:
    """
    Reads an SDF file containing one or more molecules, and for each molecule:
    - Extracts the SMILES string
    - Extracts all user-defined properties

    Args:
        sdf_file (str | Path): Path to the input SDF.

    Returns:
        List[Dict[str, Any]]: One dict per molecule with keys "smiles" and "properties".
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

        output.append({"smiles": smiles_str, "properties": properties})

    return output

def load_reference_molecule(
    ref_arg: str,
    random_seed: int = 42,
) -> Chem.Mol:
    """
    Loads the reference molecule (either from a file path or a SMILES string).
    If it has no 3D conformer, embed and optimize it. Returns a 3D-ref Mol.
    
    Args:
        ref_arg (str): Path to .sdf/.mol/.smi or a SMILES string.
        random_seed (int): Seed for ETKDG embedding if needed.

    Returns:
        Chem.Mol: 3D‐embedded reference Mol with AtomMapNum set to each atom’s index+1.
    """
    if os.path.exists(ref_arg):
        ext = os.path.splitext(ref_arg)[-1].lower()
        if ext in [".sdf", ".mol"]:
            ref = Chem.MolFromMolFile(ref_arg, removeHs=False)
            if ref is None:
                raise ValueError(f"Could not read a molecule from '{ref_arg}'.")
        elif ext in [".smi", ".smiles"]:
            with open(ref_arg, "r") as f:
                first_line = f.read().splitlines()[0].strip().split()[0]
            ref = Chem.MolFromSmiles(first_line)
            if ref is None:
                raise ValueError(f"Could not parse SMILES from '{ref_arg}'.")
        else:
            raise ValueError(f"Reference file '{ref_arg}' has unsupported extension '{ext}'.")
    else:
        ref = Chem.MolFromSmiles(ref_arg)
        if ref is None:
            raise ValueError(f"Could not interpret '{ref_arg}' as a SMILES string.")

    # Ensure the reference has one 3D conformer:
    if ref.GetNumConformers() == 0 or not ref.GetConformer().Is3D():
        ref = generate_3d(ref, random_seed=random_seed)

    ref = Chem.RemoveHs(ref)
    for idx_atom, atom in enumerate(ref.GetAtoms()):
        atom.SetAtomMapNum(idx_atom + 1)
    return ref


#@beartype
def read_sdf_properties(sdf_file: str | Path) -> dict:
    """Reads all user-defined properties from an SDF file (single molecule) and returns them as a dictionary.

    Args:
        sdf_file (str | Path): Path to the SDF.

    Returns:
        Dict[str, str]: Property name → value.
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


#@beartype
def get_properties_in_sdf_file(sdf_file: str | Path) -> list:
    """Returns a list of all user-defined properties in an SDF file

     Args:
        sdf_file (str | Path): Path to the SDF.

    Returns:
        List[str]: Unique property names.
    """
    from rdkit import Chem

    # Load molecules from the SDF file
    supplier = Chem.SDMolSupplier(str(sdf_file), sanitize=False)

    properties = []

    for _, mol in enumerate(supplier):
        if mol is None:
            continue  # Skip invalid molecules

        for prop_name in mol.GetPropNames():
            properties.append(prop_name)

    return list(set(properties))


#@beartype
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


#@beartype
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


#@beartype
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


#@beartype
def smiles_list_to_base64_png_list(
    smiles_list: list[str],
    *,
    size: tuple[int, int] = (300, 100),
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
    from rdkit.Chem import rdDepictor
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


#@beartype
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


#@beartype
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


#@beartype
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


#@beartype
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


#@beartype
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


#@beartype
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


#@beartype
def show_molecules_in_sdf_file(sdf_file: str | Path):
    """show molecules in an SDF file in a Jupyter notebook using molstar"""

    from deeporigin_molstar import JupyterViewer, MoleculeViewer

    molecule_viewer = MoleculeViewer(
        data=str(sdf_file),
        format="sdf",
    )
    html_content = molecule_viewer.render_ligand()
    JupyterViewer.visualize(html_content)


def mcs(mols: list[Chem.Mol], *, timeout: int = 10) -> Chem.Mol:
    """
    Generate the Maximum Common Substructure (MCS) for molecules

    Args:
        mols (List[Chem.Mol]]): Molecules to compare.
        timeout (int): Seconds before MCS search times out.

    Returns:
        Chem.Mol: Mol built from the MCS SMARTS string.
    """

    prepped = [preprocess_mol(m) for m in mols]

    params = rdFMCS.MCSParameters()
    params.AtomTyper = rdFMCS.AtomCompare.CompareElements
    params.BondTyper = rdFMCS.BondCompare.CompareOrder
    params.BondCompareParameters.RingMatchesRingOnly = False
    params.BondCompareParameters.CompleteRingsOnly = False
    params.Timeout = timeout
    params.Verbose = False

    result = rdFMCS.FindMCS(prepped, parameters=params)

    if result.canceled:
        raise RuntimeError("MCS computation timed out or failed.")

    return Chem.MolFromSmarts(result.smartsString)


def preprocess_mol(mol: Chem.Mol) -> Chem.Mol:
    """
    Preprocess a molecule for MCS

    Args:
        mol (Chem.Mol): RDKit molecule

    Returns:
        Chem.Mol: Preprocessed molecule
    """
    mol = Chem.RemoveHs(mol)
    Chem.SanitizeMol(mol)
    return mol

def generate_3d(mol: Chem.Mol, random_seed: int = 42) -> Chem.Mol:
    """
    Given an RDKit Mol (2D or 3D), generate a single 3D conformer via ETKDG + UFF.
    Returns a new Mol that has explicit Hs and a 3D conformer.
    
    Args:
        mol (Chem.Mol): Input 2D/3D Mol.
        random_seed (int): Seed for ETKDG embedding.

    Returns:
        Chem.Mol: 3D‐embedded Mol (explicit Hs stripped).
    """
    m = Chem.Mol(mol)
    m = Chem.AddHs(m)
    params = AllChem.ETKDGv3()
    params.randomSeed = random_seed
    success = AllChem.EmbedMolecule(m, params)
    if success != 0:
        # Retry once with a slightly different seed
        params.randomSeed = random_seed + 1
        success = AllChem.EmbedMolecule(m, params)
        if success != 0:
            raise RuntimeError("3D embedding failed for molecule.")
    AllChem.UFFOptimizeMolecule(m)
    m = Chem.RemoveHs(m)
    return m

def _extract_mcs_fragment_from_reference(
    ref_mol: Chem.Mol,
    target: Chem.Mol,
    mcs_smarts: str,
    debug: bool = False
) -> Chem.Mol:
    """
    Given a 3D reference molecule (`ref_mol`), a target molecule (`target`),
    and an MCS SMARTS string, do the following:
      1) Find the atom‐index mapping for the MCS between ref and target.
      2) Transfer each matched atom's 3D coordinates from ref onto target.
      3) Print out before/after positions for each matched atom in the target (if debug=True).
      4) Extract the MCS “core” from the reference (ReplaceSidechains + DeleteSubstructs), and return that core.

    Note: This function *mutates* `target`’s conformer (creating one if necessary).
    Args:
        ref_mol (Chem.Mol): 3D reference Mol.
        target (Chem.Mol): Target Mol to transfer onto.
        mcs_smarts (str): SMARTS pattern of the MCS.
        debug (bool): Print debug messages.

    Returns:
        Chem.Mol: Core fragment Mol extracted from reference.
    """
    if debug:
        print(f"  [_extract_mcs_fragment] SMARTS: {mcs_smarts}")

    # 1) Build an RDKit Mol from the MCS SMARTS
    mcs_pat = Chem.MolFromSmarts(mcs_smarts)
    if mcs_pat is None:
        raise ValueError("Invalid MCS SMARTS string.")

    # 2) Remove Hs on both ref and target for a clean substructure match
    ref_no_h = Chem.RemoveHs(ref_mol)
    ref_no_h.UpdatePropertyCache()
    target_no_h = Chem.RemoveHs(target)
    target_no_h.UpdatePropertyCache()

    # 3) Find the matching atom indices in ref_no_h and target_no_h
    match_ref = ref_no_h.GetSubstructMatch(mcs_pat)
    if not match_ref:
        raise ValueError(
            f"MCS pattern DID NOT match the reference molecule.\n"
            f" SMARTS: {mcs_smarts}\n"
            f" Reference (no Hs): {Chem.MolToSmiles(ref_no_h)}"
        )
    match_tgt = target_no_h.GetSubstructMatch(mcs_pat)
    if not match_tgt:
        raise ValueError(
            f"MCS pattern DID NOT match the target molecule.\n"
            f" SMARTS: {mcs_smarts}\n"
            f" Target (no Hs): {Chem.MolToSmiles(target_no_h)}"
        )

    # 4) Grab the 3D conformer from ref_mol and ensure it is 3D
    if ref_mol.GetNumConformers() == 0 or not ref_mol.GetConformer(0).Is3D():
        raise RuntimeError("Reference mol does not have a valid 3D conformer.")
    ref_conf = ref_mol.GetConformer(0)

    # 5) Ensure `target` has at least one conformer (even if it’s 2D)
    if target.GetNumConformers() == 0:
        new_conf_tgt = Chem.Conformer(target.GetNumAtoms())
        target.AddConformer(new_conf_tgt, assignId=True)
    tgt_conf = target.GetConformer(0)

    # 6) Transfer coordinates, printing debug if requested
    if debug:
        print("    Transferring coordinates for matched atoms:")
    for j, ref_idx in enumerate(match_ref):
        tgt_idx = match_tgt[j]
        ref_pos = ref_conf.GetAtomPosition(ref_idx)
        if debug:
            try:
                old_tgt_pos = tgt_conf.GetAtomPosition(tgt_idx)
                print(f"      Target atom {tgt_idx} before:  "
                      f"({old_tgt_pos.x:.3f}, {old_tgt_pos.y:.3f}, {old_tgt_pos.z:.3f})")
            except Exception:
                print(f"      Target atom {tgt_idx} before:  (no valid conformer)")
        tgt_conf.SetAtomPosition(tgt_idx, ref_pos)
        if debug:
            new_tgt_pos = tgt_conf.GetAtomPosition(tgt_idx)
            print(f"      Target atom {tgt_idx} after:   "
                  f"({new_tgt_pos.x:.3f}, {new_tgt_pos.y:.3f}, {new_tgt_pos.z:.3f}) "
                  f"(copied from Ref atom {ref_idx})")

    # 7) Extract the MCS “core” from ref_no_h via ReplaceSidechains → DeleteSubstructs
    replaced = AllChem.ReplaceSidechains(target, mcs_pat)
    if replaced is None:
        raise RuntimeError(
            "ReplaceSidechains returned None even though SMARTS matched the reference."
        )
    dummy = Chem.MolFromSmiles("*")
    core = AllChem.DeleteSubstructs(replaced, dummy)
    if core is None:
        raise RuntimeError(
            "DeleteSubstructs returned None. Either there was no '*' or something failed."
        )
    core.UpdatePropertyCache()

    if debug:
        print(f"    Core extracted (SMILES): {Chem.MolToSmiles(core)}")

    return core

def constrained_minimize(
    mol: Chem.Mol,
    atom_indices: List[int],
    tol: float = 0.05,
    force_const: float = 200.0,
    max_cycles: int = 10,
    max_iter_per_cycle: int = 1000
) -> None:
    """
    Given a 3D‐conformer on `mol`, add MMFF position‐constraints on
    each atom in `atom_indices` and run a constrained MMFF minimize.
    In other words, it holds specified atoms fixed, relax the rest.
    Modifies `mol` in place.
    
    Args:
        mol (Chem.Mol): 3D Mol to minimize.
        atom_indices (List[int]): 0‑based atoms to constrain.
        tol (float): Tolerance for positional constraints.
        force_const (float): Force constant for constraints.
        max_cycles (int): How many minimize cycles.
        max_iter_per_cycle (int): Iterations per cycle.
    """
    # build MMFF force‐field for this molecule
    mmff_props = rdForceFieldHelpers.MMFFGetMoleculeProperties(mol)
    ff = rdForceFieldHelpers.MMFFGetMoleculeForceField(mol, mmff_props)
    # add a “spring” constraint on each core atom
    for idx in atom_indices:
        ff.MMFFAddPositionConstraint(idx, tol, force_const)
    # Run minimization in chunks to ensure convergence within the specified cycles.
    cycles = max_cycles
    while ff.Minimize(maxIts=max_iter_per_cycle) and cycles > 0:
        cycles -= 1

def pairwise_align_and_annotate(
    ref_mol: Chem.Mol,
    targets: List[Chem.Mol],
    minimize: bool,
    debug: bool = False
) -> Tuple[Chem.Mol, List[Chem.Mol], List[Chem.Mol]]:
    """
    For each target in `targets`:
      1) Compute pairwise MCS vs. ref_mol.
      2) Attempt to ConstrainedEmbed (copy ref coords → target matched atoms, then UFF).
         If ConstrainedEmbed fails, catch the exception, print it (if debug=True),
         and add it to the failed list (with an Alignment_Error property).
         Then perform a normal 3D embed for that target so it still ends up with coordinates.
      3) Annotate every target (successful or not) with these SD‐properties:
         - Alignment_Success (True/False)
         - Total_Atoms
         - Target_Murcko_Atoms
         - Matched_Atoms
         - Unmatched_Atoms
         - Matched_In_Murcko
         - (if failure) Alignment_Error

    Args:
        ref_mol (Chem.Mol): Reference molecule with a valid 3D conformer.
        targets (List[Chem.Mol]]): Target molecules to align.
        minimize (bool): Whether to apply constrained optimization.
        debug (bool): Whether to print detailed debug messages.

    Returns:
        Tuple[
            Chem.Mol,          # ref_annotated: reference with per‐target annotations
            List[Chem.Mol],    # annotated_targets: aligned & annotated targets
            List[Chem.Mol]     # failed_targets: targets that could not be aligned
        ]
    """
    ref_annotated = Chem.Mol(ref_mol)  # copy for adding per‐target props
    annotated_targets: List[Chem.Mol] = []
    failed_targets: List[Chem.Mol] = []

    # Precompute reference Murcko scaffold atom indices
    ref_no_h = Chem.RemoveHs(ref_mol)
    ref_no_h.UpdatePropertyCache()
    try:
        ref_murcko = MurckoScaffold.GetScaffoldForMol(ref_no_h)
        ref_murcko.UpdatePropertyCache()
        ref_scaffold_match = ref_no_h.GetSubstructMatch(ref_murcko)
        ref_scaffold_indices = set(ref_scaffold_match)
    except Exception:
        ref_scaffold_indices = set()

    for idx, tgt_original in enumerate(targets):
        # Always remove Hs from the original target
        tgt_original = Chem.RemoveHs(tgt_original)

        if debug:
            print(f"\n--- Processing target idx {idx} ---")
            print(f"  Reference SMILES: {Chem.MolToSmiles(ref_mol)}")
            print(f"  Target    SMILES: {Chem.MolToSmiles(tgt_original)}")

        # 1) Compute pairwise MCS
        pair_mcs = rdFMCS.FindMCS(
            [ref_mol, tgt_original],
            threshold=0.9,
            ringMatchesRingOnly=True,
            atomCompare=rdFMCS.AtomCompare.CompareAny,
            bondCompare=rdFMCS.BondCompare.CompareAny,
            timeout=300,
        )
        mcs_smarts_i = pair_mcs.smartsString

        # 2) Annotate reference: which atoms match this target’s MCS?
        patt_ref = Chem.MolFromSmarts(mcs_smarts_i)
        match_ref = ref_mol.GetSubstructMatch(patt_ref)
        if match_ref:
            mapping_pairs_ref = [f"{j}:{match_ref[j]}" for j in range(len(match_ref))]
            ref_annotated.SetProp(f"MatchedAtomMap_target{idx}", ",".join(mapping_pairs_ref))
            ref_annotated.SetProp(f"MCS_target{idx}_SMARTS", mcs_smarts_i)
        else:
            ref_annotated.SetProp(f"MatchedAtomMap_target{idx}", "")
            ref_annotated.SetProp(f"MCS_target{idx}_SMARTS", "")

        # 3) Prepare target and find its match_tgt
        tgt = Chem.Mol(tgt_original)  # fresh copy
        patt_tgt = Chem.MolFromSmarts(mcs_smarts_i)
        match_tgt = tgt.GetSubstructMatch(patt_tgt)

        # Count atoms for properties
        total_atoms = tgt.GetNumAtoms()
        # Build target’s Murcko scaffold indices
        tgt_no_h = Chem.RemoveHs(tgt)
        tgt_no_h.UpdatePropertyCache()
        try:
            tgt_murcko = MurckoScaffold.GetScaffoldForMol(tgt_no_h)
            tgt_murcko.UpdatePropertyCache()
            tgt_scaffold_match = tgt_no_h.GetSubstructMatch(tgt_murcko)
            tgt_scaffold_indices = set(tgt_scaffold_match)
        except Exception:
            tgt_scaffold_indices = set()

        num_tgt_murcko_atoms = len(tgt_scaffold_indices)
        num_matched = len(match_tgt)
        num_unmatched = total_atoms - num_matched

        # Count how many matched atoms lie in both scaffolds
        num_matched_in_murcko = 0
        for j, ref_idx in enumerate(match_ref):
            tgt_idx = match_tgt[j]
            if (ref_idx in ref_scaffold_indices) and (tgt_idx in tgt_scaffold_indices):
                num_matched_in_murcko += 1

        num_unmatched_murcko_atoms = num_tgt_murcko_atoms - num_matched_in_murcko
        murcko_coverage = num_matched_in_murcko * 2 / (len(ref_scaffold_indices) + len(tgt_scaffold_indices))

        # 4) Attempt alignment
        alignment_success = True
        error_msg = ""
        if not match_tgt:
            # No MCS → normal 3D embed
            tgt_3d = generate_3d(tgt, random_seed=42 + idx)
            alignment_success = False
            error_msg = "No pairwise MCS match."
            if debug:
                print(f"  [Alignment] No pairwise MCS match for target {idx}")
        else:
            try:
                core_frag = _extract_mcs_fragment_from_reference(
                    ref_mol, tgt, mcs_smarts_i, debug=debug
                )
                tgt_3d = AllChem.ConstrainedEmbed(
                    tgt,
                    core_frag
                )
                if debug:
                    print(f"  [Alignment] Success for target {idx}")
            except Exception as ex:
                alignment_success = False
                error_msg = str(ex)
                if debug:
                    print(f"  [Alignment] Failed for target {idx}: {error_msg}")
                failed_copy = Chem.Mol(tgt)  # 2D-only tgt
                failed_copy.SetProp("Alignment_Error", error_msg)
                # Then give it a normal 3D so it still ends up with coordinates
                #failed_copy = generate_3d(failed_copy, random_seed=42 + idx)
                failed_targets.append(failed_copy)

        if alignment_success:
            # 5) Annotate the target (tgt_3d) with SD properties
            tgt_3d.SetProp("Alignment_Success", str(alignment_success))
            tgt_3d.SetProp("Total_Atoms", str(total_atoms))
            tgt_3d.SetProp("Target_Murcko_Atoms", str(num_tgt_murcko_atoms))
            tgt_3d.SetProp("Matched_Atoms", str(num_matched))
            tgt_3d.SetProp("Unmatched_Atoms", str(num_unmatched))
            tgt_3d.SetProp("Matched_In_Murcko", str(num_matched_in_murcko))
            tgt_3d.SetProp("UnMatched_In_Murcko", str(num_unmatched_murcko_atoms))
            tgt_3d.SetProp("Murcko_coverage", str(murcko_coverage))

            # 6) Add atom‐mapping numbers so you can visualize which atoms matched
            mapping_pairs_tgt: List[str] = []
            for pat_idx, tgt_idx in enumerate(match_tgt):
                ref_idx = match_ref[pat_idx] + 1
                atom = tgt_3d.GetAtomWithIdx(tgt_idx)
                # now using the raw ref_idx
                atom.SetAtomMapNum(ref_idx)
                mapping_pairs_tgt.append(f"{ref_idx}:{tgt_idx}")
            if mapping_pairs_tgt:
                tgt_3d.SetProp("MatchedAtomMap", ",".join(mapping_pairs_tgt))

            # 7) Add Name if necessary
            if not tgt_3d.HasProp("_Name") or not tgt_3d.GetProp("_Name").strip():
                tgt_3d.SetProp("_Name", f"Target_{idx}")

            # 8) Minimized if necessary
            if minimize:
                constrained_minimize(tgt_3d, list(match_tgt))
        else:
            tgt_3d.SetProp("Alignment_Error", error_msg)

        annotated_targets.append(tgt_3d)

    return ref_annotated, annotated_targets, failed_targets

def consolidate_reference_annotations(
    ref_annotated: Chem.Mol
) -> Chem.Mol:
    """
    Gather all per-target MCS annotations on the reference into two consolidated properties:
      - MatchedAtomMap_All:  "targetIdx:smartsIdx:refIdx,smartsIdx:refIdx|targetIdx2:…"
      - MCS_All_SMARTS:       "targetIdx1:SMARTS1|targetIdx2:SMARTS2|…"
    Returns the same `ref_annotated` molecule with two new props added.
    """
    all_map_entries: List[str] = []
    all_smarts_entries: List[str] = []

    for prop_name in ref_annotated.GetPropNames():
        if prop_name.startswith("MatchedAtomMap_target"):
            target_idx = prop_name.replace("MatchedAtomMap_target", "")
            mapping = ref_annotated.GetProp(prop_name)
            all_map_entries.append(f"{target_idx}:{mapping}")
        elif prop_name.startswith("MCS_target") and prop_name.endswith("_SMARTS"):
            target_idx = prop_name.replace("MCS_target", "").replace("_SMARTS", "")
            smarts = ref_annotated.GetProp(prop_name)
            all_smarts_entries.append(f"{target_idx}:{smarts}")

    if all_map_entries:
        ref_annotated.SetProp("MatchedAtomMap_All", "|".join(all_map_entries))
    if all_smarts_entries:
        ref_annotated.SetProp("MCS_All_SMARTS", "|".join(all_smarts_entries))

    return ref_annotated

def run_pairwise_alignment(
    input_mols: List[Chem.Mol],
    reference: Chem.Mol,
    minimize: bool,
    debug: bool
) -> list[Chem.Mol]:
    """
    Align a list of molecules to a reference via pairwise MCS alignment.

    Args:
        input_mols (List[Chem.Mol]]): Ligand molecules to align.
        reference (Chem.Mol): Reference molecule (with AtomMapNum set and 3D coords).
        minimize (bool): If True, apply constrained MMFF minimization to each ligand.
        debug (bool): If True, print debug info and mapped‐SMILES.

    Returns:
        List[Chem.Mol]: The list of successfully aligned (and annotated) molecules.
    """

    #add mapping
    for atom in reference.GetAtoms():
        atom.SetAtomMapNum(atom.GetIdx() + 1)

    print("Performing pairwise MCS, alignment, and annotation...")
    ref_annotated, annotated_targets, failed_records = pairwise_align_and_annotate(
        reference,
        input_mols,
        minimize,
        debug=debug
    )

    # 4) Consolidate per-target properties on the reference
    #ref_consolidated = consolidate_reference_annotations(ref_annotated)

    if debug:
        for i, tgt in enumerate(annotated_targets):
            print(f"Target {i} SMILES with maps:")
            print(Chem.MolToSmiles(tgt, isomericSmiles=True))
        print(f"Reference SMILES with maps:")
        print(Chem.MolToSmiles(reference, isomericSmiles=True))

    return annotated_targets

def align(
    *,
    mols: list[Chem.Mol],
    reference: Chem.Mol,
    energy: float = 5,
    minimize: bool,
    debug: bool
) -> list[list[dict]]:
    """
    Aligns a set of molecules to a reference and returns MCS atom constraints.

    Args:
        mols (List[Chem.Mol]]): Ligands to align.
        reference (Chem.Mol): Reference molecule (3D, AtomMapNum set).
        mcs_mol (Chem.Mol): SMARTS‐derived MCS for alignment.
        energy (float): Constraint energy weight.
        minimize (bool): Whether to apply constrained MMFF minimization.
        debug (bool): Print debug info (mapped SMILES, etc.).

    Returns:
        List[List[Dict[str, Any]]]]: For each ligand, a list of dicts:
            - "index" (int): target atom’s AtomMapNum
            - "coordinates" (List[float]): [x,y,z] from reference
            - "energy" (float): energy weight
    """

    #align Molecules
    aligned_mols = run_pairwise_alignment(
        mols,
        reference,
        minimize,
        debug
    )

    # Generate constraints
    all_constraints = []

    for i, mol in enumerate(aligned_mols):
        # Build constraints: atom index + 1 (1-based), and position from ref
        constraints = [
            {
                "index": atom.GetAtomMapNum(),  
                "coordinates": list(mol.GetConformer(0).GetAtomPosition(idx)),
                "energy": energy,
            }
            for idx, atom in enumerate(mol.GetAtoms())
        ]

        all_constraints.append(constraints)

    return all_constraints