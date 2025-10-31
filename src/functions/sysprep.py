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
    retain_waters: bool = False,
    add_H_atoms: bool = True,
    protonate_protein: bool = True,
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
        "add_H_atoms": add_H_atoms,
        "protonate_protein": protonate_protein,
        "retain_waters": retain_waters,
        "padding": padding,
        "use_cache": use_cache,
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
    # response = tools_api.run_function(
    #     key="deeporigin.system-prep",
    #     version="0.3.2",
    #     function_execution_params_schema_dto=body,
    # )

    # mock it for now
    # TODO -- remove once function is deployed
    response = {
        "status": "success",
        "message": "System preparation completed successfully",
        "protein_path": "entities/proteins/db4aa32e2e8ffa976a60004a8361b86427a2e5653a6623bb60b7913445902549.pdb",
        "ligand_path": "entities/ligands/bac7b4d01c1a7ab102b1c9955a1839730a5099b08eba93807e12f6ab22adfb67.sdf",
        "output_files": [
            "function-runs/system-prep/586a55e3205e7b807c49b32d085700e55e4c159b672f32f7a276ad1b1cc1e9d1/bsm_system.xml",
            "function-runs/system-prep/586a55e3205e7b807c49b32d085700e55e4c159b672f32f7a276ad1b1cc1e9d1/system.pdb",
            "function-runs/system-prep/586a55e3205e7b807c49b32d085700e55e4c159b672f32f7a276ad1b1cc1e9d1/solvation.xml",
        ],
        "function_key": "deeporigin.system-prep",
        "function_version": "0.3.2",
        "parameters": {
            "add_H_atoms": False,
            "protonate_protein": False,
            "retain_waters": True,
            "padding": 1.5,
        },
    }

    prepared_system_pdb_path = [
        file for file in response["output_files"] if file.endswith(".pdb")
    ][0]

    file_api.download_file(
        remote_path=prepared_system_pdb_path,
        local_path=output_pdb_path,
    )

    return output_pdb_path
