"""This module contains the function to run loop modelling on a protein identified by a PDB ID."""

import hashlib
import os

from beartype import beartype
import requests

# URL = "http://loop-modelling.default.jobs.edge.deeporigin.io/model_loops"
URL = "http://localhost:8080/model_loops"
CACHE_DIR = os.path.expanduser("~/.deeporigin/model_loops")


@beartype
def model_loops(
    *,
    pdb_id: str,
    use_cache: bool = True,
) -> dict:
    """
    Run system preparation on a protein-ligand complex.

    Args:
        pdb_id: PDB ID of the protein to model loops for.

    Returns:
        Path to the output PDB file if successful, or raises RuntimeError if the server fails.
    """

    # Create a hash of the input parameters for caching
    hash_input = f"{pdb_id}"
    cache_key = hashlib.md5(hash_input.encode()).hexdigest()
    cache_path = os.path.join(CACHE_DIR, cache_key)
    output_json_path = os.path.join(cache_path, "data.json")
    output_pdb_path = os.path.join(cache_path, "data.pdb")

    # Create cache directory if it doesn't exist
    os.makedirs(cache_path, exist_ok=True)

    # Check if cached results exist
    if use_cache and os.path.exists(output_json_path):
        import json

        with open(output_json_path, "r") as f:
            return json.load(f)

    # If no cached results, proceed with server call
    payload = {
        "pdb_id": pdb_id,
    }

    try:
        response = requests.post(URL, json=payload, stream=True)
    except Exception as e:
        raise RuntimeError(f"Failed to connect to server: {e}") from None

    if response.status_code == 200:
        import json

        response_json = response.json()
        # Write pdb_contents to output_pdb_path
        pdb_contents = response_json.pop("pdb_contents", None)
        if pdb_contents is not None:
            with open(output_pdb_path, "w") as pdb_file:
                pdb_file.write(pdb_contents)

        response_json["pdb_file"] = output_pdb_path
        # Write the modified JSON (without pdb_contents) to output_json_path
        with open(output_json_path, "w") as f:
            json.dump(response_json, f, indent=2)
        return response_json

    # If the server request fails, raise an error
    raise RuntimeError(
        f"Server returned status code {response.status_code}: {response.text}"
    ) from None
