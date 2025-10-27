"""this module implements a low level function to find pockets in a protein determined by a PDB file"""

import os

from deeporigin.drug_discovery.structures import Protein
from deeporigin.utils.core import hash_dict

URL = "http://pocketfinder.default.jobs.edge.deeporigin.io/find_pockets"
CACHE_DIR = os.path.expanduser("~/.deeporigin/pocket-finder")


# @beartype
def find_pockets(
    *,
    protein: Protein,
    pocket_count: int = 5,
    pocket_min_size: int = 30,
    use_cache: bool = True,
) -> str | None:
    """Find protein binding pockets in a PDB structure and save the results.

    This function sends a PDB file to a remote server for pocket detection and
    saves the returned results to a cache directory.

    Args:
        protein (Protein): protein to find pockets in
        pocket_count (int, optional): Maximum number of pockets to detect. Defaults to 5.
        pocket_min_size (int, optional): Minimum size of pockets to consider. Defaults to 30.

    Returns:
        str | None: Path to the cache directory if successful, None if the request failed.
    """

    if pocket_count < 1:
        raise ValueError("pocket_count must be at least 1") from None
    if pocket_min_size < 1:
        raise ValueError("pocket_min_size must be at least 1") from None

    # Prepare the request payload
    payload = {
        "protein_path": protein._remote_path,
        "pocket_count": pocket_count,
        "pocket_min_size": pocket_min_size,
    }

    cache_key = hash_dict(payload)
    cache_path = os.path.join(CACHE_DIR, cache_key)

    # Check if cached results exist
    if use_cache and os.path.exists(cache_path):
        return cache_path

    protein.upload()
    os.makedirs(cache_path, exist_ok=True)

    from deeporigin.platform import file_api, tools_api

    body = {"params": payload, "clusterId": tools_api.get_default_cluster_id()}

    # Send the request to the server
    response = tools_api.run_function(
        key="deeporigin.pocketfinder",
        version="0.2.0",
        function_execution_params_schema_dto=body,
    )

    for file in response.files:
        file_api.download_file(
            remote_path=file,
            local_path=os.path.join(cache_path, file.split("/")[-1]),
        )

    return cache_path
