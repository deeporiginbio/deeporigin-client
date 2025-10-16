"""This module contains the function to run system preparation on a protein-ligand complex."""

import os

from beartype import beartype

from deeporigin.drug_discovery.structures import Ligand, Protein
from deeporigin.utils.core import hash_dict

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
        "protein_path": protein._remote_path,
        "ligand_path": ligand._remote_path,
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

    protein.upload()
    ligand.upload()

    # If no cached results, proceed with server call
    from deeporigin.platform import file_api, tools_api

    body = {"params": payload, "clusterId": tools_api.get_default_cluster_id()}

    # Send the request to the server
    response = tools_api.run_function(
        key="deeporigin.system-prep",
        version="0.2.0",
        body=body,
    )

    file_api.download_file(
        remote_path=response.system_pdb_path,
        local_path=output_pdb_path,
    )

    return output_pdb_path
