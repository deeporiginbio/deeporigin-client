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
from typing import Optional, Tuple

import requests
from beartype import beartype

from deeporigin.drug_discovery.structures import Pocket, Protein
from deeporigin.exceptions import DeepOriginException

URL = "http://docking.default.jobs.edge.deeporigin.io/dock"
CACHE_DIR = os.path.expanduser("~/.deeporigin/docking")


@beartype
def dock(
    protein: Protein,
    smiles_string: str,
    box_size: Tuple[float, float, float] = (20.0, 20.0, 20.0),
    pocket_center: Optional[Tuple[int, int, int]] = None,
    pocket: Optional[Pocket] = None,
):
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

    # Create hash of inputs
    hasher = hashlib.sha256()

    # Hash protein file contents
    with open(protein.file, "rb") as f:
        hasher.update(f.read())

    # Hash pocket file if provided
    if pocket is not None and pocket.file is not None:
        with open(pocket.file, "rb") as f:
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
        with open(protein.file, "rb") as f:
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
        with open(sdf_file, "w") as file:
            for solution in response[0]["solutions"]:
                file.write(solution["output_sdf_content"])

    from deeporigin_molstar import DockingViewer, JupyterViewer

    docking_viewer = DockingViewer()
    html_content = docking_viewer.render_with_seperate_crystal(
        protein_data=str(protein.file),
        protein_format="pdb",
        ligands_data=[sdf_file],
        ligand_format="sdf",
    )

    JupyterViewer.visualize(html_content)
