"""
ligand_alignment.py

This module provides utilities for aligning ligands to a reference pose using a Maximum Common Substructure (MCS) approach.
It includes functions to extract atom coordinates for the MCS in a reference molecule, map test ligands to these coordinates,
and generate constraints for use in constrained docking workflows.


Functions:
    - getRefMCSCoordMap: Extracts coordinates of MCS atoms from a reference molecule.
    - createTestMaps: Maps a test ligand to the MCS and generates coordinate constraints.
    - align_ligands_to_reference: Main entry point to align a set of ligands to a reference pose using an MCS SMARTS string.
"""

from beartype import beartype
from rdkit import Chem


def getRefMCSCoordMap(referencePose, MCSMol):
    """
    Extracts the 3D coordinates of atoms in the reference molecule that match the Maximum Common Substructure (MCS).

    Args:
        referencePose (rdkit.Chem.Mol): Reference molecule with 3D coordinates (e.g., from an SD file).
        MCSMol (rdkit.Chem.Mol): RDKit molecule representing the MCS (from SMARTS).

    Returns:
        list[list[float]]: List of [x, y, z] coordinates for each atom in the reference that matches the MCS.

    Raises:
        ValueError: If the MCS SMARTS does not match the reference molecule.
    """

    # Get MCS atom indices in reference molecule
    refMatch = referencePose.GetSubstructMatch(MCSMol)

    if not refMatch:
        raise ValueError("MCS SMARTS did not match the reference molecule.")

    # Get the x/y/z coordinates of the corresponding atoms
    refPosePositions = [
        list(referencePose.GetConformer().GetAtomPosition(atom_idx))
        for i, atom_idx in enumerate(refMatch)
    ]

    return refPosePositions


def createTestMaps(refPosePositions, SMILES, MCSMol):
    """
    Maps a test ligand (by SMILES) to the MCS and generates coordinate constraints for each matching atom.

    Args:
        refPosePositions (list[list[float]]): Coordinates of MCS atoms in the reference pose.
        SMILES (str): SMILES string of the test ligand to be aligned.
        MCSMol (rdkit.Chem.Mol): RDKit molecule representing the MCS (from SMARTS).

    Returns:
        list[dict]: List of constraints, each a dict with keys 'index' (1-based atom index),
                    'coordinates' ([x, y, z]), and 'energy' (restraint strength).
    """
    mol = Chem.MolFromSmiles(SMILES)
    matchingIndexes = mol.GetSubstructMatch(MCSMol)

    constraints = [
        {"index": x + 1, "coordinates": refPosePositions[i], "energy": 5}
        for i, x in enumerate(matchingIndexes)
    ]

    return constraints


@beartype
def align_ligands_to_reference(
    *,
    smiles_strings: list[str],
    mcs_smarts_string: str,
    reference_pose_filepath: str,
):
    """
    Aligns a set of ligands to a reference pose using a Maximum Common Substructure (MCS) SMARTS string.

    For each test ligand, finds the MCS match, and generates coordinate constraints to align the ligand to the reference pose.

    Args:
        smiles_strings (list[str]): List of SMILES strings for ligands to be aligned (must contain the MCS).
        mcs_smarts_string (str): SMARTS string representing the MCS present in all ligands.
        reference_pose_filepath (str): Filepath to the reference pose (SDF file with 3D coordinates).

    Returns:
        list[list[dict]]: For each test ligand, a list of constraint dicts describing atom indices (1-based) and target coordinates.
    """

    # Create MCSMol
    MCSMol = Chem.MolFromSmarts(mcs_smarts_string)

    # Get the reference pose
    referencePose = Chem.SDMolSupplier(reference_pose_filepath)[0]

    # Get a dict that describes {atom_index : [x/y/z]} for MCS in reference molecule
    refMap = getRefMCSCoordMap(referencePose, MCSMol)

    # Array of dicts that map atom_index to x/y/z coordinates they should be constrained to (refMap)
    testMaps = [createTestMaps(refMap, x, MCSMol) for x in smiles_strings]

    return testMaps
