"""this module contains low-level functions to interact
with the data API. functions here simply wrap API endpoints."""

import os
from dataclasses import dataclass
from typing import Literal, Optional, Union

import requests
from beartype import beartype
from deeporigin import cache_do_api_tokens, get_do_api_tokens
from deeporigin.config import get_value
from deeporigin.exceptions import DeepOriginException
from deeporigin.utils import _nucleus_url

# types of rows
RowType = Literal["row", "workspace", "database"]


@dataclass
class DeepOriginClient:
    api_url = _nucleus_url()
    org_id = get_value()["organization_id"]

    api_access_token, api_refresh_token = get_do_api_tokens()
    cache_do_api_tokens(api_access_token, api_refresh_token)

    headers = {
        "accept": "application/json",
        "authorization": f"Bearer {api_access_token}",
        "content-type": "application/json",
        "x-org-id": org_id,
    }

    @beartype
    def invoke(
        self,
        endpoint: str,
        data: dict,
    ) -> Union[dict, list]:
        """core call to API endpoint"""

        response = requests.post(
            f"{self.api_url}{endpoint}",
            headers=self.headers,
            json=data,
        )

        return _check_response(response)


# default client
CLIENT = DeepOriginClient()


@beartype
def list_rows(
    *,
    parent_id: Optional[str] = None,
    row_type: Optional[RowType] = None,
    parent_is_root: Optional[bool] = None,
    client: DeepOriginClient = CLIENT,
) -> list[dict]:
    """low level API that wraps the ListRows endpoint"""

    filters = []

    if parent_is_root is not None:
        filters.append(dict(parent=dict(isRoot=parent_is_root)))

    if parent_id:
        filters.append(dict(parent=dict(id=parent_id)))

    if row_type:
        filters.append(dict(rowType=row_type))

    data = dict(filters=filters)
    return client.invoke("ListRows", data)


@beartype
def describe_database_stats(
    database_id: str,
    *,
    client: DeepOriginClient = CLIENT,
):
    return client.invoke("DescribeDatabaseStats", dict(databaseId=database_id))


@beartype
def list_mentions(
    query: str,
    *,
    client: DeepOriginClient = CLIENT,
):
    return client.invoke("ListMentions", dict(query=query))


@beartype
def list_row_back_references(row_id: str, *, client: DeepOriginClient = CLIENT):
    return client.invoke("ListRowBackReferences", dict(rowId=row_id))


@beartype
def create_file_download_url(
    file_id: str,
    *,
    client: DeepOriginClient = CLIENT,
) -> dict:
    """low-level API call to CreateFileDownloadUrl"""
    return client.invoke("CreateFileDownloadUrl", dict(fileId=file_id))


@beartype
def describe_file(
    file_id: str,
    *,
    client: DeepOriginClient = CLIENT,
) -> dict:
    """low-level API call to DescribeFile"""
    return client.invoke("DescribeFile", dict(fileId=file_id))


@beartype
def describe_row(
    row_id: str,
    *,
    fields: bool = False,
    client: DeepOriginClient = CLIENT,
) -> dict:
    """low-level API that wraps the DescribeRow endpoint."""

    data = dict(rowId=row_id, fields=fields)
    return client.invoke("DescribeRow", data)


@beartype
def list_database_rows(
    row_id: str,
    *,
    client: DeepOriginClient = CLIENT,
) -> list[dict]:
    """low level API that wraps the ListDatabaseRows endpoint"""

    data = dict(databaseRowId=row_id)
    return client.invoke("ListDatabaseRows", data)


@beartype
def download_file(
    file_id: str, destination: str, *, client: DeepOriginClient = CLIENT
) -> None:
    """download the file to the destination folder"""

    if not os.path.isdir(destination):
        raise DeepOriginException(f"{destination} should be a path to a folder.")

    file_name = describe_file(file_id, client=client)["name"]

    url = create_file_download_url(file_id, client=client)["downloadUrl"]

    save_path = os.path.join(destination, file_name)

    response = requests.get(url)
    if response.status_code == 200:
        with open(save_path, "wb") as file:
            file.write(response.content)
    else:
        raise DeepOriginException(f"Failed to download file {file_id}")


@beartype
def convert_id_format(
    *,
    hids: Optional[Union[list[str], set[str]]] = None,
    ids: Optional[Union[list[str], set[str]]] = None,
    client: DeepOriginClient = CLIENT,
) -> list[dict]:
    """convert a list of HIDs to IDs or vice versa"""

    if hids is None and ids is None:
        raise DeepOriginException(
            "Either `hids` or `ids` should be non-None and a list of strings"
        )

    conversions = []

    if hids is not None:
        for hid in hids:
            conversions.append(dict(hid=hid))

    if ids is not None:
        for sid in ids:
            conversions.append(dict(id=sid))

    data = dict(conversions=conversions)

    return client.invoke("ConvertIdFormat", data)


@beartype
def _check_response(response: requests.models.Response) -> Union[dict, list]:
    """utility function to check responses"""

    if response.status_code == 404:
        raise DeepOriginException("[Error 404] The requested resource was not found.")

    response.raise_for_status()
    response = response.json()

    if "error" in response:
        raise DeepOriginException(response["error"])

    if "data" in response:
        return response["data"]
    else:
        raise KeyError("`data` not in response")
