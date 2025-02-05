"""module that contains some utility functions for chemistry"""

from typing import Optional

from beartype import beartype
from rdkit import Chem


@beartype
def smiles_to_sdf(smiles: str, sdf_path: str) -> None:
    """convert a SMILES string to a SDF file"""

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
def sdf_to_smiles(sdf_file: str) -> Optional[str]:
    """Extracts the first molecule's SMILES string from an SDF file using RDKit.

    Args:
        sdf_file (str): Path to the SDF file.

    Returns:
        Optional[str]: SMILES string if extraction is successful, else None.
    """

    suppl = Chem.SDMolSupplier(sdf_file)
    if not suppl or suppl[0] is None:
        return None

    return Chem.MolToSmiles(suppl[0])


@beartype
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
