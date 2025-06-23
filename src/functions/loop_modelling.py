"""This module contains the function to run loop modelling on a protein identified by a PDB ID."""

import hashlib
import os

from beartype import beartype
import requests

URL = "http://localhost:8080/model_loops"
URL = "http://loop-modelling.default.jobs.edge.deeporigin.io/model_loops"
CACHE_DIR = os.path.expanduser("~/.deeporigin/model_loops")


@beartype
def model_loops(
    *,
    pdb_id: str,
) -> str:
    """
    Run system preparation on a protein-ligand complex.

    Args:
        protein_path (str | Path): Path to the protein file.
        ligand_path (str | Path): Path to the ligand file.
        padding (float, optional): Padding to add around the system. Defaults to 1.0.
        keep_waters (bool, optional): Whether to keep water molecules. Defaults to False.
        is_lig_protonated (bool, optional): Whether the ligand is already protonated. Defaults to True.
        is_protein_protonated (bool, optional): Whether the protein is already protonated. Defaults to True.

    Returns:
        Path to the output PDB file if successful, or raises RuntimeError if the server fails.
    """

    # Create a hash of the input parameters for caching
    hash_input = f"{pdb_id}"
    cache_key = hashlib.md5(hash_input.encode()).hexdigest()
    cache_path = os.path.join(CACHE_DIR, cache_key)
    output_pdb_path = os.path.join(cache_path, "loops_modelled.pdb")

    # Check if cached results exist
    if os.path.exists(output_pdb_path):
        return output_pdb_path

    # If no cached results, proceed with server call
    payload = {
        "pdb_id": pdb_id,
    }

    response = requests.post(URL, json=payload, stream=True)

    if response.status_code == 200:
        # Create cache directory if it doesn't exist
        os.makedirs(cache_path, exist_ok=True)
        # Save the result to the cache
        with open(output_pdb_path, "wb") as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)
        return output_pdb_path

    # If the server request fails, raise an error
    raise RuntimeError(
        f"Server returned status code {response.status_code}: {response.text}"
    ) from None
