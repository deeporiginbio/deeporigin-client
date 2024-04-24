"""this module contains low-level functions to interact
with the data API. functions here simply wrap API endpoints."""

import os
from typing import Optional, Union

import requests
from beartype import beartype
from deeporigin.exceptions import DeepOriginException
from deeporigin.managed_data.client import Client, DeepOriginClient
from deeporigin.managed_data.schema import RowType


def _get_default_client(client):
    if client is None:
        client = DeepOriginClient()  # pragma: no cover
        client.authenticate()  # pragma: no cover
    return client


@beartype
def list_rows(
    *,
    parent_id: Optional[str] = None,
    row_type: Optional[RowType] = None,
    parent_is_root: Optional[bool] = None,
    client: Optional[Client] = None,
) -> list[dict]:
    """low level API that wraps the ListRows endpoint"""

    client = _get_default_client(client)
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
def list_files(
    *,
    assigned_row_ids: Optional[list[str]] = None,
    is_unassigned: Optional[bool] = None,
    client: Optional[Client] = None,
):
    client = _get_default_client(client)

    filters = []

    if is_unassigned is not None:
        filters.append(dict(isUnassigned=is_unassigned))

    if assigned_row_ids:
        filters.append(dict(assignedRowIds=assigned_row_ids))

    return client.invoke("ListFiles", data=dict(filters=filters))


@beartype
def describe_database_stats(
    database_id: str,
    *,
    client: Optional[Client] = None,
):
    client = _get_default_client(client)

    return client.invoke("DescribeDatabaseStats", dict(databaseId=database_id))


@beartype
def list_mentions(
    query: str,
    *,
    client: Optional[Client] = None,
):
    client = _get_default_client(client)

    return client.invoke("ListMentions", dict(query=query))


@beartype
def list_row_back_references(
    row_id: str,
    *,
    client: Optional[Client] = None,
):
    client = _get_default_client(client)

    return client.invoke("ListRowBackReferences", dict(rowId=row_id))


@beartype
def create_file_download_url(
    file_id: str,
    *,
    client: Optional[Client] = None,
) -> dict:
    """low-level API call to CreateFileDownloadUrl"""

    client = _get_default_client(client)

    return client.invoke("CreateFileDownloadUrl", dict(fileId=file_id))


@beartype
def describe_file(
    file_id: str,
    *,
    client: Optional[Client] = None,
) -> dict:
    """low-level API call to DescribeFile"""

    client = _get_default_client(client)

    return client.invoke("DescribeFile", dict(fileId=file_id))


@beartype
def describe_row(
    row_id: str,
    *,
    fields: bool = False,
    client: Optional[Client] = None,
) -> dict:
    """low-level API that wraps the DescribeRow endpoint."""

    client = _get_default_client(client)

    data = dict(rowId=row_id, fields=fields)

    return client.invoke("DescribeRow", data)


@beartype
def list_database_rows(
    row_id: str,
    *,
    client: Optional[Client] = None,
) -> list[dict]:
    """low level API that wraps the ListDatabaseRows endpoint"""

    client = _get_default_client(client)

    data = dict(databaseRowId=row_id)
    return client.invoke("ListDatabaseRows", data)


@beartype
def download_file(
    file_id: str,
    *,
    destination: str = os.getcwd(),
    client: Optional[Client] = None,
) -> None:
    """download the file to the destination folder"""

    client = _get_default_client(client)

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
    client: Optional[Client] = None,
) -> list[dict]:
    """convert a list of HIDs to IDs or vice versa"""

    if hids is None and ids is None:
        raise DeepOriginException(
            "Either `hids` or `ids` should be non-None and a list of strings"
        )

    client = _get_default_client(client)

    conversions = []

    if hids is not None:
        for hid in hids:
            conversions.append(dict(hid=hid))

    if ids is not None:
        for sid in ids:
            conversions.append(dict(id=sid))

    data = dict(conversions=conversions)

    return client.invoke("ConvertIdFormat", data)
