import base64
import requests
from beartype import beartype
from typing import Tuple, Optional
from deeporigin.drug_discovery.structures import Pocket, Protein
from deeporigin.exceptions import DeepOriginException

URL = "http://docking.default.jobs.edge.deeporigin.io/dock"


@beartype
def dock(
    protein: Protein,
    smiles_list: list[str],
    box_size: Tuple[float, float, float] = (20.0, 20.0, 20.0),
    pocket_center: Optional[Tuple[int, int, int]] = None,
    pocket: Optional[Pocket] = None,
):
    """
    Run molecular docking using the DeepOrigin API.

    Args:
        protein (Protein): Protein object representing the target protein
        smiles_list (list[str]): List of SMILES strings for ligands
        box_size (Tuple[float, float, float]): Size of the docking box (x, y, z)
        pocket_center (Tuple[int, int, int]): Center coordinates of the docking pocket (x, y, z)

    Returns:
        dict: API response
    """

    if pocket is not None:
        pocket_center = pocket.get_center().tolist()

    if pocket_center is None:
        raise DeepOriginException("Pocket center is required")

    # Read and encode the protein file
    with open(protein.file, "rb") as f:
        encoded_protein = base64.b64encode(f.read()).decode("utf-8")

    # Prepare the request payload
    payload = {
        "functionId": "docking",
        "params": {
            "protein": encoded_protein,
            "smiles_list": smiles_list,
            "box_size": list(box_size),
            "pocket_center": list(pocket_center),
        },
    }

    # Make the API request
    response = requests.post(
        URL,
        json=payload["params"],
        headers={"Content-Type": "application/json"},
    )

    # Raise an exception for bad status codes
    response.raise_for_status()

    response = response.json()

    # return response

    sdf_file = "temp.sdf"

    with open(sdf_file, "w") as file:
        for solution in response[0]["solutions"]:
            file.write(solution["output_sdf_content"])

    from deeporigin_molstar import DockingViewer, JupyterViewer

    docking_viewer = DockingViewer()
    html_content = docking_viewer.render_with_seperate_crystal(
        protein_data=str(protein.file),
        protein_format="pdb",
        ligands_data=[sdf_file],
        ligand_format="sdf",
    )

    JupyterViewer.visualize(html_content)
