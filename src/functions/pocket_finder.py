"""this module implements a low level function to find pockets in a protein determined by a PDB file"""

import base64
import os
import zipfile

from beartype import beartype
import requests

from deeporigin.drug_discovery.structures import Protein
from deeporigin.utils.core import hash_dict

URL = "http://pocketfinder.default.jobs.edge.deeporigin.io/find_pockets"
CACHE_DIR = os.path.expanduser("~/.deeporigin/pocket-finder")


@beartype
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

    # Prepare the request payload
    payload = {
        "protein": protein.to_base64(),
        "pocket_count": pocket_count,
        "pocket_min_size": pocket_min_size,
    }

    cache_key = hash_dict(payload)
    cache_path = os.path.join(CACHE_DIR, cache_key)

    # Check if cached results exist
    if use_cache and os.path.exists(cache_path):
        return cache_path

    # Send the request to the server
    response = requests.post(URL, json=payload)

    # Check if the request was successful
    if response.status_code == 200:
        # Create cache directory if it doesn't exist
        os.makedirs(cache_path, exist_ok=True)

        # Decode base64 content and save as zip file
        zip_path = os.path.join(cache_path, "results.zip")
        with open(zip_path, "wb") as f:
            f.write(base64.b64decode(response.text))

        # Extract the zip file to the cache directory
        with zipfile.ZipFile(zip_path, "r") as zip_ref:
            zip_ref.extractall(cache_path)

        # Remove the zip file
        os.remove(zip_path)

        return cache_path
    print(f"Error: Server returned status code {response.status_code}")
    print(response.text)
    return None
