"""this module contains low-level functions to interact
with the data API. functions here simply wrap API endpoints."""

import os
from typing import Optional, Union

import requests
from beartype import beartype
from deeporigin import cache_do_api_tokens, get_do_api_tokens
from deeporigin.config import get_value
from deeporigin.exceptions import DeepOriginException

API_URL = get_value()["nucleus_api_endpoint"]
ORG_ID = get_value()["organization_id"]


@beartype
def download_file(file_id: str, destination: str) -> None:
    """download the file to the destination folder"""

    if not os.path.isdir(destination):
        raise DeepOriginException(f"{destination} should be a path to a folder.")

    file_name = describe_file(file_id)["name"]

    url = create_file_download_url(file_id)["downloadUrl"]

    save_path = os.path.join(destination, file_name)

    response = requests.get(url)
    if response.status_code == 200:
        with open(save_path, "wb") as file:
            file.write(response.content)
    else:
        raise DeepOriginException(f"Failed to download file {file_id}")


@beartype
def create_file_download_url(file_id: str) -> dict:
    """low-level API call to CreateFileDownloadUrl"""
    return _invoke("CreateFileDownloadUrl", dict(fileId=file_id))


@beartype
def describe_file(file_id: str) -> dict:
    """low-level API call to DescribeFile"""
    return _invoke("DescribeFile", dict(fileId=file_id))


@beartype
def describe_row(row_id: str, *, fields: bool = False) -> dict:
    """low-level API that wraps the DescribeRow endpoint."""

    data = dict(rowId=row_id, fields=fields)
    return _invoke("DescribeRow", data)


@beartype
def list_rows(row_id: str) -> list[dict]:
    """low level API that wraps the ListRows endpoint"""

    data = dict(filters=[dict(parent=dict(id=row_id))])
    return _invoke("ListRows", data)


@beartype
def list_database_rows(row_id: str) -> list[dict]:
    """low level API that wraps the ListDatabaseRows endpoint"""

    data = dict(databaseRowId=row_id)
    return _invoke("ListDatabaseRows", data)


@beartype
def convert_id_format(
    *,
    hids: Optional[list[str]] = None,
    ids: Optional[list[str]] = None,
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

    return _invoke("ConvertIdFormat", data)


@beartype
def _invoke(endpoint: str, data: dict) -> Union[dict, list]:
    """core call to API endpoint"""

    response = requests.post(
        f"{API_URL}{endpoint}",
        headers=_generate_headers(),
        json=data,
    )

    return _check_response(response)


@beartype
def _generate_headers() -> dict:
    """generate headers used for requests to the API"""

    api_access_token, api_refresh_token = get_do_api_tokens()
    cache_do_api_tokens(api_access_token, api_refresh_token)

    headers = {
        "accept": "application/json",
        "authorization": f"Bearer {api_access_token}",
        "content-type": "application/json",
        "x-org-id": ORG_ID,
    }

    return headers


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
