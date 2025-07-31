"""This module implements a low level function to perform molecular docking using the Deep Origin API.

The main function `dock()` takes a Protein object, a list of ligand SMILES strings, and docking box parameters
to perform docking calculations. The docking box can be specified either by providing explicit coordinates for the
pocket center, or by passing a Pocket object which contains the pocket center information.

The module interfaces with the Deep Origin docking service to perform the actual docking calculations remotely.
"""

import base64
import hashlib
import json
import os
from pathlib import Path
from typing import Optional
import zipfile

from beartype import beartype
from rdkit import Chem
import requests

from deeporigin.drug_discovery.structures import Ligand, Pocket, Protein
from deeporigin.exceptions import DeepOriginException


@beartype
def dock(
    *,
    protein: Protein,
    smiles_string: Optional[str] = None,
    ligand: Optional[Ligand] = None,
    box_size: tuple[float, float, float] = (20.0, 20.0, 20.0),
    pocket_center: Optional[tuple[int, int, int]] = None,
    pocket: Optional[Pocket] = None,
    use_cache: bool = True,
) -> str:
    """
    Run molecular docking using the DeepOrigin API.

    Args:
        protein (Protein): Protein object representing the target protein
        smiles_list (list[str]): List of SMILES strings for ligands
        box_size (Tuple[float, float, float]): Size of the docking box (x, y, z)
        pocket_center (Tuple[int, int, int]): Center coordinates of the docking pocket (x, y, z)

    Returns:
        dict: API response
    """

    URL = "http://docking.default.jobs.edge.deeporigin.io/dock"
    CACHE_DIR = os.path.expanduser("~/.deeporigin/docking")

    if pocket is not None:
        pocket_center = pocket.get_center().tolist()

    if pocket_center is None:
        raise DeepOriginException("Pocket center is required") from None

    if ligand is not None:
        smiles_string = ligand.smiles

    if smiles_string is None:
        raise DeepOriginException(
            "Either smiles_string or ligand must be provided"
        ) from None

    # Create hash of inputs
    hasher = hashlib.sha256()

    # Hash protein file contents
    protein_file = protein._dump_state()
    with open(protein_file, "rb") as f:
        hasher.update(f.read())

    # Hash pocket file if provided
    if pocket is not None and pocket.file_path is not None:
        with open(pocket.file_path, "rb") as f:
            hasher.update(f.read())

    # Hash other inputs
    hasher.update(
        json.dumps(
            {
                "smiles": smiles_string,
                "box_size": list(box_size),
                "pocket_center": list(pocket_center),
            }
        ).encode()
    )

    cache_hash = hasher.hexdigest()
    sdf_file = str(Path(CACHE_DIR) / f"{cache_hash}.sdf")

    # Check if cached result exists
    if os.path.exists(sdf_file) and use_cache:
        return sdf_file
    else:
        # Read and encode the protein file
        with open(protein_file, "rb") as f:
            encoded_protein = base64.b64encode(f.read()).decode("utf-8")

        # Prepare the request payload
        payload = {
            "functionId": "docking",
            "params": {
                "protein": encoded_protein,
                "smiles_list": [smiles_string],
                "box_size": list(box_size),
                "pocket_center": list(pocket_center),
            },
        }

        # Make the API request
        response = requests.post(
            URL,
            json=payload["params"],
            headers={"Content-Type": "application/json"},
        )

        # Raise an exception for bad status codes
        response.raise_for_status()

        response = response.json()

        # Write SDF file to cache
        Path(sdf_file).parent.mkdir(parents=True, exist_ok=True)
        with open(sdf_file, "w") as file:
            for solution in response[0]["solutions"]:
                file.write(solution["output_sdf_content"])

    return sdf_file


def constrained_dock(
    *,
    protein: Protein,
    reference_ligand: Ligand,
    ligand: Ligand,
    mcs: Chem.Mol,
    pocket: Optional[Pocket] = None,
    box_size: tuple[float, float, float] = (20.0, 20.0, 20.0),
    pocket_center: Optional[tuple[int, int, int]] = None,
    use_cache: bool = True,
) -> list[str]:
    """Perform constrained molecular docking using a reference ligand and MCS alignment.

    This function performs molecular docking with constraints based on a reference ligand
    and a maximum common substructure (MCS). The ligand is aligned to the reference ligand
    using the MCS, and docking is performed with these alignment constraints applied.

    The function uses a remote docking service and implements caching to avoid redundant
    computations. Results are cached based on a hash of all input parameters.

    Args:
        protein: The protein structure to dock against.
        reference_ligand: The reference ligand used for alignment constraints.
        ligand: The ligand to be docked.
        mcs: The maximum common substructure (RDKit Mol object) between reference and target ligands.
        pocket: Optional pocket object. If provided, its center will be used as pocket_center.
        box_size: Size of the docking box in Angstroms (x, y, z). Defaults to (20.0, 20.0, 20.0).
        pocket_center: Optional center coordinates for the docking box. If None and pocket is provided,
                     uses the pocket center. If both are None, raises an error.
        use_cache: Whether to use cached results if available. Defaults to True.

    Returns:
        list[str]: List of file paths to the docking result files (typically SDF files).

    Raises:
        ValueError: If neither pocket_center nor pocket is provided.
        requests.RequestException: If the docking service request fails.
        zipfile.BadZipFile: If the response from the service is not a valid zip file.

    Note:
        The function creates a cache directory at ~/.deeporigin/constrained_docking/ and
        stores results based on a SHA256 hash of all input parameters. This allows for
        efficient reuse of previous docking results.
    """
    URL = "http://constrained-docking.default.jobs.edge.deeporigin.io/dock"
    CACHE_DIR = os.path.expanduser("~/.deeporigin/constrained_docking")

    from deeporigin.drug_discovery import chemistry

    constraints = chemistry.align(
        mols=[ligand.mol.m],
        reference=reference_ligand.mol.m,
        mcs_mol=mcs,
    )[0]

    # get pocket center
    if pocket_center is None:
        pocket_center = pocket.get_center().tolist()

    # Create hash of inputs
    hasher = hashlib.sha256()

    # Hash protein file contents
    protein_file = protein._dump_state()
    with open(protein_file, "rb") as f:
        hasher.update(f.read())

    # Hash reference ligand
    reference_ligand_b64 = reference_ligand.to_base64()
    hasher.update(reference_ligand_b64.encode())

    # Hash ligand
    ligand_b64 = ligand.to_base64()
    hasher.update(ligand_b64.encode())

    # Hash other inputs
    hasher.update(
        json.dumps(
            {
                "constraints": constraints,
                "box_size": list(box_size),
                "pocket_center": list(pocket_center),
                "top_criteria": "energy_score",
            }
        ).encode()
    )

    cache_hash = hasher.hexdigest()
    zip_file = str(Path(CACHE_DIR) / f"{cache_hash}.zip")
    extract_dir = str(Path(CACHE_DIR) / cache_hash)

    # Check if cached result exists
    if os.path.exists(extract_dir) and use_cache:
        # Return paths to extracted files
        extracted_files = []
        for file_path in Path(extract_dir).glob("*"):
            if file_path.is_file():
                extracted_files.append(str(file_path))
        return extracted_files
    else:
        payload = {
            "box_size": box_size,
            "constraints": constraints,
            "ligands": "/data/ligand.sdf",
            "protein": {
                "file_path": "/data/protein.pdb",
                "pocket_center": pocket_center,
            },
            "top_criteria": "energy_score",
        }

        # Inject the base64-encoded strings into the payload
        payload["protein_b64"] = protein.to_base64()
        payload["ligand_b64"] = ligand.to_base64()

        # Send the POST request
        response = requests.post(
            URL,
            headers={"Content-Type": "application/json"},
            json=payload,
        )

        # Raise an exception for bad status codes
        response.raise_for_status()

        # Write zip file to cache temporarily
        Path(zip_file).parent.mkdir(parents=True, exist_ok=True)
        with open(zip_file, "wb") as f:
            f.write(response.content)

        # Extract zip file
        Path(extract_dir).mkdir(parents=True, exist_ok=True)
        with zipfile.ZipFile(zip_file, "r") as zip_ref:
            zip_ref.extractall(extract_dir)

        # Delete the zip file
        os.remove(zip_file)

        # Return paths to extracted files
        extracted_files = []
        for file_path in Path(extract_dir).glob("*"):
            if file_path.is_file():
                extracted_files.append(str(file_path))
        return extracted_files
