"""This module contains the function to run system preparation on a protein-ligand complex."""

import os

from beartype import beartype
import requests

from deeporigin.drug_discovery.structures import Ligand, Protein
from deeporigin.utils.core import hash_dict

URL = "http://sysprep.default.jobs.edge.deeporigin.io/sysprep"
CACHE_DIR = os.path.expanduser("~/.deeporigin/sysprep")


@beartype
def run_sysprep(
    *,
    protein: Protein,
    ligand: Ligand,
    padding: float = 1.0,
    keep_waters: bool = False,
    is_lig_protonated: bool = True,
    is_protein_protonated: bool = True,
    use_cache: bool = True,
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

    payload = {
        "protein": protein.to_base64(),
        "ligand": ligand.to_base64(),
        "is_lig_protonated": is_lig_protonated,
        "is_protein_protonated": is_protein_protonated,
        "keep_waters": keep_waters,
        "padding": padding,
    }

    # Create a hash of the input parameters for caching

    cache_key = hash_dict(payload)
    cache_path = os.path.join(CACHE_DIR, cache_key)
    output_pdb_path = os.path.join(cache_path, "complex.pdb")

    # Check if cached results exist
    if os.path.exists(output_pdb_path) and use_cache:
        return output_pdb_path

    # If no cached results, proceed with server call

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
