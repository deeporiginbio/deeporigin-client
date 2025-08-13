"""This module implements a low level function to perform molecular docking using the Deep Origin API.

The main function `dock()` takes a Protein object, a list of ligand SMILES strings, and docking box parameters
to perform docking calculations. The docking box can be specified either by providing explicit coordinates for the
pocket center, or by passing a Pocket object which contains the pocket center information.

The module interfaces with the Deep Origin docking service to perform the actual docking calculations remotely.
"""

import concurrent.futures
import os
from pathlib import Path
import time
from typing import Optional
import zipfile

from beartype import beartype
import requests

from deeporigin.drug_discovery.structures import Ligand, LigandSet, Pocket, Protein
from deeporigin.exceptions import DeepOriginException
from deeporigin.utils.core import hash_dict


def _extract_cached_files(extract_dir: str) -> list[str]:
    """Extract and return paths to cached files."""

    extracted_files = []
    for file_path in Path(extract_dir).glob("*"):
        if file_path.is_file():
            extracted_files.append(str(file_path))
    return extracted_files


def _get_pocket_center(
    pocket: Optional[Pocket],
    pocket_center: Optional[tuple[int, int, int]],
) -> list:
    """Extract pocket center coordinates."""
    if pocket_center is None:
        return pocket.get_center().tolist()
    return list(pocket_center)


@beartype
def dock(
    *,
    protein: Protein,
    smiles_string: Optional[str] = None,
    ligand: Optional[Ligand] = None,
    box_size: tuple[float, float, float] = (20.0, 20.0, 20.0),
    pocket_center: Optional[tuple[int, int, int]] = None,
    pocket: Optional[Pocket] = None,
    use_cache: bool = True,
) -> str:
    """
    Run molecular docking using the DeepOrigin API.

    Args:
        protein (Protein): Protein object representing the target protein
        smiles_string (Optional[str]): SMILES string for the ligand to dock
        ligand (Optional[Ligand]): Ligand object to dock
        box_size (tuple[float, float, float]): Size of the docking box (x, y, z)
        pocket_center (Optional[tuple[int, int, int]]): Center coordinates of the docking pocket (x, y, z)
        pocket (Optional[Pocket]): Pocket object defining the docking region
        use_cache (bool): Whether to use cached results. Defaults to True

    Returns:
        str: path to the SDF file containing the docking results
    """

    URL = "http://docking.default.jobs.edge.deeporigin.io/dock"
    CACHE_DIR = os.path.expanduser("~/.deeporigin/docking")

    if pocket is not None or pocket_center is not None:
        pocket_center = _get_pocket_center(pocket, pocket_center)
    else:
        raise DeepOriginException("Pocket center is required") from None

    if ligand is not None:
        smiles_string = ligand.smiles

    if smiles_string is None:
        raise DeepOriginException(
            "Either smiles_string or ligand must be provided"
        ) from None

    # Prepare the request payload
    payload = {
        "protein": protein.to_base64(),
        "smiles_list": [smiles_string],
        "box_size": list(box_size),
        "pocket_center": list(pocket_center),
    }

    cache_hash = hash_dict(payload)
    sdf_file = str(Path(CACHE_DIR) / f"{cache_hash}.sdf")

    # Check if cached result exists
    if os.path.exists(sdf_file) and use_cache:
        return sdf_file
    else:
        # Make the API request
        response = requests.post(
            URL,
            json=payload,
            headers={"Content-Type": "application/json"},
        )

        # Raise an exception for bad status codes
        response.raise_for_status()

        response = response.json()

        # Write SDF file to cache
        Path(sdf_file).parent.mkdir(parents=True, exist_ok=True)
        with open(sdf_file, "w") as file:
            for solution in response[0]["solutions"]:
                file.write(solution["output_sdf_content"])

    return sdf_file


def constrained_dock(
    *,
    protein: Protein,
    ligand: Ligand,
    constraints: list[dict],
    pocket: Optional[Pocket] = None,
    box_size: tuple[float, float, float] = (20.0, 20.0, 20.0),
    pocket_center: Optional[tuple[int, int, int]] = None,
    use_cache: bool = True,
) -> list[str]:
    """Perform constrained molecular docking using a reference ligand constraints.

    This function performs molecular docking with constraints based on a reference ligand
    and a maximum common substructure (MCS). The ligand is aligned to the reference ligand
    using the MCS, and docking is performed with these alignment constraints applied.

    Args:
        protein (Protein): The protein structure to dock against.
        ligand (Ligand): The ligand to be docked.
        constraints (list[dict]): List of constraints for the docking. Generate this using `align.compute_constraints`.
        pocket (Optional[Pocket]): Optional pocket object. If provided, its center will be used as pocket_center.
        box_size (tuple[float, float, float]): Size of the docking box in Angstroms (x, y, z). Defaults to (20.0, 20.0, 20.0).
        pocket_center (Optional[tuple[int, int, int]]): Optional center coordinates for the docking box. If None and pocket is provided,
                     uses the pocket center. If both are None, raises an error.
        use_cache (bool): Whether to use cached results if available. Defaults to True.

    Returns:
        list[str]: List of file paths to the docking result files (typically SDF files).


    Note:
        The function creates a cache directory at ~/.deeporigin/constrained_docking/ and
        stores results based on a SHA256 hash of all input parameters. This allows for
        efficient reuse of previous docking results.
    """
    URL = "http://constrained-docking.default.jobs.edge.deeporigin.io/dock"
    CACHE_DIR = os.path.expanduser("~/.deeporigin/constrained_docking")

    if pocket is None and pocket_center is None:
        raise DeepOriginException(
            "Either pocket or pocket_center must be provided"
        ) from None

    pocket_center_list = _get_pocket_center(pocket, pocket_center)

    payload = {
        "box_size": box_size,
        "constraints": constraints,
        "protein": {"pocket_center": pocket_center_list},
        "top_criteria": "score",
        "protein_b64": protein.to_base64(),
        "ligand_b64": ligand.to_base64(),
    }

    cache_hash = hash_dict(payload)
    zip_file = str(Path(CACHE_DIR) / f"{cache_hash}.zip")
    extract_dir = str(Path(CACHE_DIR) / cache_hash)

    if os.path.exists(extract_dir) and use_cache:
        return _extract_cached_files(extract_dir)

    # Send the POST request
    response = requests.post(
        URL,
        headers={"Content-Type": "application/json"},
        json=payload,
    )
    response.raise_for_status()

    # Write and extract zip file
    Path(zip_file).parent.mkdir(parents=True, exist_ok=True)
    with open(zip_file, "wb") as f:
        f.write(response.content)

    Path(extract_dir).mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(zip_file, "r") as zip_ref:
        zip_ref.extractall(extract_dir)

    os.remove(zip_file)
    return _extract_cached_files(extract_dir)


def _parallel_dock(
    *,
    protein: Protein,
    ligands: LigandSet,
    pocket: Pocket,
    batch_size: int = 10,
    max_retries: int = 3,
    sleep_between_batches: float = 0.1,
    use_cache: bool = True,
) -> dict[str, object]:
    """
    Run docking in parallel batches with retries for failures.

    Args:
        protein: protein argument to be passed to dock().
        ligands: A ligand set to be docked
        pocket: A pocket in the protein to dock into
        batch_size: number of ligands to process in parallel.
        max_retries: max times to retry failed dockings.
        sleep_between_batches: delay between batches, in seconds.
        use_cache: whether to use cached results if available. Defaults to True.

    Returns:
        Dict with:
            - output_paths: list of output paths or None for failures
            - total_failures: total number of failures encountered (including retries)
            - permanent_failures: indices of ligands that failed even after retries
            - elapsed_time: time in seconds for entire process
            - durations: list of time spent docking each ligand (or None if permanently failed)
    """
    total = len(ligands)
    output_paths = [None] * total
    retries_left = [max_retries] * total
    total_failures = 0
    durations = [None] * total  # Track time spent per ligand

    def dock_ligand_timed(idx: int) -> Optional[tuple]:
        nonlocal total_failures
        try:
            start = time.time()
            result = dock(
                ligand=ligands[idx],
                pocket=pocket,
                protein=protein,
                use_cache=use_cache,
            )
            end = time.time()
            return (idx, result, end - start)
        except Exception as e:
            print(f"[Error] Ligand {idx} failed: {e}")
            total_failures += 1
            return (idx, None, None)

    start_time = time.time()

    while any(
        (out is None and retries > 0)
        for out, retries in zip(output_paths, retries_left, strict=False)
    ):
        to_process = [
            i
            for i, (out, retries) in enumerate(
                zip(output_paths, retries_left, strict=False)
            )
            if out is None and retries > 0
        ]
        print(f"Processing batch of {min(batch_size, len(to_process))} ligands...")

        for i in range(0, len(to_process), batch_size):
            batch = to_process[i : i + batch_size]
            with concurrent.futures.ThreadPoolExecutor() as executor:
                results = list(executor.map(dock_ligand_timed, batch))
            for idx, result, duration in results:
                if result is not None:
                    output_paths[idx] = result
                    durations[idx] = duration
                else:
                    retries_left[idx] -= 1
            time.sleep(sleep_between_batches)

    elapsed_time = time.time() - start_time
    permanent_failures = [i for i, out in enumerate(output_paths) if out is None]

    return {
        "output_paths": output_paths,
        "total_failures": total_failures,
        "permanent_failures": permanent_failures,
        "elapsed_time": elapsed_time,
        "durations": durations,
    }
