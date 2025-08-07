"""This module implements a low level function to perform molecular property predictions
using the DeepOrigin API.
"""

import json
import os
from pathlib import Path

from beartype import beartype
import requests

from deeporigin.utils.core import hash_dict

URL = "http://molprops.default.jobs.edge.deeporigin.io/properties"
CACHE_DIR = os.path.expanduser("~/.deeporigin/molprops")

# Ensure cache directory exists
os.makedirs(CACHE_DIR, exist_ok=True)


@beartype
def molprops(
    smiles_string: str,
    *,
    use_cache: bool = True,
) -> dict:
    """
    Run molecular property prediction using the DeepOrigin API.

    Args:
        smiles_string (str): SMILES string for the molecule
        use_cache (bool): Whether to use the cache

    Returns:
        str: Path to the cached SDF file containing the results
    """

    payload = {
        "smiles": smiles_string,
    }

    # Create hash of inputs
    cache_hash = hash_dict(payload)
    response_file = str(Path(CACHE_DIR) / f"{cache_hash}.json")

    # Check if cached result exists
    if os.path.exists(response_file) and use_cache:
        # Read cached response
        with open(response_file, "r") as file:
            response = json.load(file)

    else:
        # Prepare the request payload

        # Make the API request
        response = requests.post(
            URL,
            json=payload,
            headers={"Content-Type": "application/json"},
        )

        # Raise an exception for bad status codes
        response.raise_for_status()

        response = response.json()

        # Write JSON response to cache
        with open(response_file, "w") as file:
            json.dump(response, file)

    return response
