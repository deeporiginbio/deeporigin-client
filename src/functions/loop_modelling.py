"""This module contains the function to run loop modelling on a protein identified by a PDB ID."""

import os

from beartype import beartype
import requests

from deeporigin.utils.core import hash_dict

URL = "http://loop-modelling.default.jobs.edge.deeporigin.io/model_loops"
CACHE_DIR = os.path.expanduser("~/.deeporigin/model_loops")


@beartype
def model_loops(
    *,
    pdb_id: str,
    use_cache: bool = True,
) -> str:
    """
    Run system preparation on a protein-ligand complex.

    Args:
        pdb_id: PDB ID of the protein to model loops for.

    Returns:
        Path to the output PDB file if successful, or raises RuntimeError if the server fails.
    """

    # If no cached results, proceed with server call
    payload = {
        "pdb_id": pdb_id,
    }

    cache_hash = hash_dict(payload)
    cache_path = os.path.join(CACHE_DIR, cache_hash)
    output_pdb_path = os.path.join(cache_path, "loops_modelled.pdb")

    # Check if cached results exist
    if use_cache and os.path.exists(output_pdb_path):
        return output_pdb_path

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
