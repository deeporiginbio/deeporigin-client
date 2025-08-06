"""this module contains functions for aligning RDKit molecules."""

import copy

import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem, rdFMCS


def randomize_mol_pose(mol: Chem.Mol, seed: int = None) -> Chem.Mol:
    """
    Takes a 3D RDKit molecule and returns a copy with a randomly rotated and translated conformation.

    Args:
        mol: RDKit molecule to randomize
        seed: Random seed for reproducibility (default: None for random seed)
    """
    mol = copy.deepcopy(mol)
    conf = mol.GetConformer()

    # Random rotation matrix
    rng = np.random.default_rng(seed)
    theta = rng.uniform(0, 2 * np.pi)
    phi = rng.uniform(0, 2 * np.pi)
    z = rng.uniform(0, 2 * np.pi)

    # Rotation matrices around the x, y, and z axes
    Rx = np.array(
        [
            [1, 0, 0],
            [0, np.cos(theta), -np.sin(theta)],
            [0, np.sin(theta), np.cos(theta)],
        ]
    )
    Ry = np.array(
        [[np.cos(phi), 0, np.sin(phi)], [0, 1, 0], [-np.sin(phi), 0, np.cos(phi)]]
    )
    Rz = np.array([[np.cos(z), -np.sin(z), 0], [np.sin(z), np.cos(z), 0], [0, 0, 1]])

    # Combined rotation
    R = Rz @ Ry @ Rx

    # Random translation vector
    translation = rng.uniform(-5, 5, size=(3,))

    # Apply to each atom
    for i in range(mol.GetNumAtoms()):
        pos = np.array(conf.GetAtomPosition(i))
        new_pos = R @ pos + translation
        conf.SetAtomPosition(i, new_pos)

    return mol


def compute_constraints(
    *,
    mols: list[Chem.Mol],
    reference: Chem.Mol,
    mcs_mol: Chem.Mol,
    energy: float = 5,
) -> list[list[dict]]:
    """
    Aligns a set of molecules to a reference and returns MCS atom constraints.

    Args:
        mols (list[Chem.Mol]): Molecules to align.
        reference (Chem.Mol): Reference molecule (with 3D coords).
        mcs_mol (Chem.Mol): MCS molecule.
        energy (float): Energy weight for constraints.

    Returns:
        list[list[dict]]: Constraints for each molecule.
    """

    ref = preprocess_mol(reference)
    mcs_match_ref = safe_substruct_match(ref, mcs_mol, "reference")
    ref_conf = ref.GetConformer()

    # Get reference atom positions for the MCS atoms
    ref_positions = [list(ref_conf.GetAtomPosition(idx)) for idx in mcs_match_ref]

    all_constraints = []

    for i, mol in enumerate(mols):
        mol_p = preprocess_mol(mol)
        match = safe_substruct_match(mol_p, mcs_mol, f"mol #{i}")

        # Align molecule to reference using MCS
        AllChem.AlignMol(
            mol_p, ref, atomMap=list(zip(match, mcs_match_ref, strict=False))
        )

        # Build constraints: atom index + 1 (1-based), and position from ref
        constraints = [
            {
                "index": atom_idx + 1,  # +1 for 1-based indexing
                "coordinates": ref_positions[i],
                "energy": energy,
            }
            for i, atom_idx in enumerate(match)
        ]

        all_constraints.append(constraints)

    return all_constraints


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


def safe_substruct_match(
    mol: Chem.Mol,
    query: Chem.Mol,
    label: str,
) -> list[int]:
    """
    Safely get a substructure match for a molecule

    Args:
        mol (Chem.Mol): RDKit molecule
        query (Chem.Mol): Query molecule
        label (str): Label for the molecule

    Returns:
        list[int]: List of atom indices that match the query
    """
    match = mol.GetSubstructMatch(query)
    if not match:
        raise ValueError(
            f"MCS does not match {label}.\n"
            f"  SMARTS: {Chem.MolToSmarts(query)}\n"
            f"  Mol SMILES: {Chem.MolToSmiles(mol)}"
        )
    return match


def mcs(mols: list[Chem.Mol], *, timeout: int = 10) -> Chem.Mol:
    """
    Generate the Maximum Common Substructure (MCS) for molecules

    Returns:
        Mol: MCS molecule constructed from the smarts string

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
