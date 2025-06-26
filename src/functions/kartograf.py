"""This module contains the function to run loop modelling on a protein identified by a PDB ID."""

import base64
import hashlib
import json
import os
from pathlib import Path
from typing import Literal

from beartype import beartype
import requests

URL = "http://ddkit-kartograf.default.jobs.edge.deeporigin.io/map_network"
URL = "http://localhost:8080/konnektor"
CACHE_DIR = os.path.expanduser("~/.deeporigin/konnektor")


@beartype
def map_network(
    *,
    sdf_file: str,
    operation: Literal["mapping", "network", "full"] = "network",
    network_type: Literal["star", "mst", "cyclic"] = "mst",
    use_cache: bool = True,
) -> dict:
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

    # Create hash of inputs
    hasher = hashlib.sha256()

    # Read and encode the ligand file
    with open(sdf_file, "rb") as f:
        encoded_file = base64.b64encode(f.read()).decode("utf-8")

    # Prepare the request payload
    payload = {
        "operation": operation,
        "file_b64": encoded_file,
        "network_type": network_type,
    }

    # Hash inputs
    hasher.update(json.dumps(payload).encode())

    cache_hash = hasher.hexdigest()
    results_json = str(Path(CACHE_DIR) / f"{cache_hash}.json")

    # Check if cached result exists
    if not os.path.exists(results_json) or not use_cache:
        # Make the API request
        response = requests.post(
            URL,
            json=payload,
            headers={"Content-Type": "application/json"},
        )

        # Raise an exception for bad status codes
        response.raise_for_status()

        response = response.json()

        # Write SDF file to cache
        Path(results_json).parent.mkdir(parents=True, exist_ok=True)
        with open(results_json, "w") as f:
            json.dump(response, f)

    else:
        with open(results_json, "r") as f:
            response = json.load(f)

    return response
