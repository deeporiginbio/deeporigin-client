"""module to run AutoDock Vina on Deep Origin"""

import functools
from pathlib import Path
from typing import Optional

import do_sdk_platform
import requests
from beartype import beartype
from deeporigin.auth import get_tokens
from deeporigin.config import get_value
from deeporigin.data_hub import api

access_token = get_tokens()["access"]
headers = {
    "authorization": f"Bearer {access_token}",
    "content-type": "application/json",
}


@functools.cache
@beartype
def _get_cluster_id() -> str:
    """gets a valid cluster ID to run tools on"""

    configuration = do_sdk_platform.configuration.Configuration(
        host="https://os.edge.deeporigin.io/api",
        access_token=access_token,
    )

    client = do_sdk_platform.ApiClient(configuration=configuration)
    client = do_sdk_platform.ClustersApi(api_client=client)

    data = client.list_clusters_controller_list_clusters(
        org_friendly_id=get_value()["organization_id"]
    )
    clusters = data.data
    cluster = [
        cluster for cluster in clusters if "us-west-2" in cluster.attributes.name
    ]
    cluster = cluster[0]
    cluster_id = cluster.id
    return cluster_id


def start_run(
    *,
    receptor: str | dict | Path = None,
    ligand: str | dict | Path = None,
    search_space: dict = None,
    docking: dict = None,
    database_id: Optional[str] = None,
    row_id: Optional[str] = None,
) -> dict:
    """starts an AutoDock Vina run"""

    if not database_id.startswith("_database"):
        data = api.resolve_ids(database_ids=[database_id])
        database_id = data[0].id

    if not row_id.startswith("_row"):
        data = api.resolve_ids(
            rows=[
                {
                    "databaseId": "vina",
                    "rowId": "vina-1",
                }
            ]
        )
        row_id = data[0].id

    json_data = {
        "toolId": "deeporigin/autodock-vina",
        "inputs": {
            "receptor": {
                "rowId": row_id,
                "columnId": "_column:NlQmxQVOF5qdlfzdgCsWw",
                "databaseId": database_id,
            },
            "ligand": {
                "rowId": row_id,
                "columnId": "_column:xoXnzZnT2bXOLBvzfcsxI",
                "databaseId": database_id,
            },
            "searchSpace": {
                "center_x": 15.190,
                "center_y": 53.903,
                "center_z": 16.917,
                "size_x": 1.1,
                "size_y": 1.1,
                "size_z": 1.1,
            },
            "docking": {
                "energy_range": 0.3,
                "exhaustiveness": 1,
                "num_modes": 9,
            },
        },
        "outputs": {
            "output_file": {
                "rowId": row_id,
                "columnId": "_column:24YqSNYeyDlpZYxrbbBdW",
                "databaseId": database_id,
            },
        },
        "clusterId": _get_cluster_id(),
    }

    response = requests.post(
        "https://os.edge.deeporigin.io/api/tools/likely-aardvark-ewo",
        headers=headers,
        json=json_data,
    )

    return response.json()


def query_run():
    pass
