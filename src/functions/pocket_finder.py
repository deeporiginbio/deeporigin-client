"""this module implements a low level function to find pockets in a protein determined by a PDB file"""

import base64
import hashlib
import os
from pathlib import Path
import zipfile

from beartype import beartype
import requests

URL = "http://pocketfinder.default.jobs.edge.deeporigin.io/find_pockets"
# URL = "http://localhost:8080/find_pockets"
CACHE_DIR = os.path.expanduser("~/.deeporigin/pocket-finder")


@beartype
def encode_pdb_file(file_path: str) -> str:
    """Read a PDB file and return its base64 encoded content."""
    with open(file_path, "rb") as f:
        return base64.b64encode(f.read()).decode("utf-8")


@beartype
def find_pockets(
    pdb_file_path: str | Path,
    pocket_count: int = 5,
    pocket_min_size: int = 30,
) -> str | None:
    """Find protein binding pockets in a PDB structure and save the results.

    This function sends a PDB file to a remote server for pocket detection and
    saves the returned results to a cache directory.

    Args:
        pdb_file_path (str | Path): Path to the PDB file to analyze.
        pocket_count (int, optional): Maximum number of pockets to detect. Defaults to 5.
        pocket_min_size (int, optional): Minimum size of pockets to consider. Defaults to 30.

    Returns:
        str | None: Path to the cache directory if successful, None if the request failed.
    """
    # Create a hash of the input parameters
    pdb_file_path = str(pdb_file_path)
    hash_input = f"{pdb_file_path}:{pocket_count}:{pocket_min_size}"
    cache_key = hashlib.md5(hash_input.encode()).hexdigest()
    cache_path = os.path.join(CACHE_DIR, cache_key)

    # Check if cached results exist
    if os.path.exists(cache_path):
        return cache_path

    # If no cached results, proceed with server call
    # Base64 encode the PDB file
    encoded_pdb = encode_pdb_file(pdb_file_path)

    # Prepare the request payload
    payload = {
        "protein": encoded_pdb,
        "pocket_count": pocket_count,
        "pocket_min_size": pocket_min_size,
    }

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
