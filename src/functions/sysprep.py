"""This module contains the function to run system preparation on a protein-ligand complex."""

import base64
import hashlib
import os
from pathlib import Path

from beartype import beartype
import requests

URL = "http://sysprep.default.jobs.edge.deeporigin.io/sysprep"
CACHE_DIR = os.path.expanduser("~/.deeporigin/sysprep")


@beartype
def sysprep(
    *,
    protein_path: str | Path,
    ligand_path: str | Path,
    padding: float = 1.0,
    keep_waters: bool = False,
    is_lig_protonated: bool = True,
    is_protein_protonated: bool = True,
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
    # Ensure paths are strings for file operations and hashing
    protein_path_str = str(protein_path)
    ligand_path_str = str(ligand_path)

    # Read and base64-encode the files
    with open(protein_path_str, "rb") as f:
        protein_b64 = base64.b64encode(f.read()).decode("utf-8")
    with open(ligand_path_str, "rb") as f:
        ligand_b64 = base64.b64encode(f.read()).decode("utf-8")

    # Create a hash of the input parameters for caching
    hash_input = f"{protein_path_str}:{ligand_path_str}:{is_lig_protonated}:{is_protein_protonated}:{keep_waters}:{padding}"
    cache_key = hashlib.md5(hash_input.encode()).hexdigest()
    cache_path = os.path.join(CACHE_DIR, cache_key)
    output_pdb_path = os.path.join(cache_path, "complex.pdb")

    # Check if cached results exist
    if os.path.exists(output_pdb_path):
        return output_pdb_path

    # If no cached results, proceed with server call
    payload = {
        "protein": protein_b64,
        "ligand": ligand_b64,
        "is_lig_protonated": is_lig_protonated,
        "is_protein_protonated": is_protein_protonated,
        "keep_waters": keep_waters,
        "padding": padding,
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
    )
