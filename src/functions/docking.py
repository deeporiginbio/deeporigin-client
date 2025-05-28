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

from beartype import beartype
import requests

from deeporigin.drug_discovery.structures import Ligand, Pocket, Protein
from deeporigin.exceptions import DeepOriginException

URL = "http://docking.default.jobs.edge.deeporigin.io/dock"
CACHE_DIR = os.path.expanduser("~/.deeporigin/docking")


@beartype
def dock(
    *,
    protein: Protein,
    smiles_string: Optional[str] = None,
    ligand: Optional[Ligand] = None,
    box_size: tuple[float, float, float] = (20.0, 20.0, 20.0),
    pocket_center: Optional[tuple[int, int, int]] = None,
    pocket: Optional[Pocket] = None,
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

    if pocket is not None:
        pocket_center = pocket.get_center().tolist()

    if pocket_center is None:
        raise DeepOriginException("Pocket center is required")

    if ligand is not None:
        smiles_string = ligand.smiles

    if smiles_string is None:
        raise DeepOriginException("Either smiles_string or ligand must be provided")

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
    if not os.path.exists(sdf_file):
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
