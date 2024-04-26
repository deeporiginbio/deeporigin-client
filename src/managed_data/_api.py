"""The `deeporigin.managed_data._api` module contains low-level functions for interacting
with Deep Origin's managed data API. The functions in this module
simply provide Pythonic interfaces to individual API endpoints."""

import os
from typing import Optional, Union

import requests
from beartype import beartype
from deeporigin.exceptions import DeepOriginException
from deeporigin.managed_data.client import Client, DeepOriginClient
from deeporigin.managed_data.schema import RowType


def _get_default_client(client: Optional[Client] = None):
    """Internal function to instantiate client

    Creates and returns an authenticated client if
    not provided with one.

    Warning: Internal function
        Do not use this function

    Args:
        client: None, or a Client


    """
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
    """Low level function that wraps the `ListRows` endpoint.

    Returns a list of rows from workspaces and databases,
    based on the parent, row type, or whether the parent is the root.

    Args:
        parent_id: ID (or human ID) or the parent.
        row_type: One of `row`, `workspace`, or `database`.
        parent_is_root: If `True` only rows that are children of the root will be returned.

    Returns:
        A list of dictionaries, where each entry corresponds to a row. Each dictionary conforms to a [ListRowsResponse][src.managed_data.schema.ListRowsResponse].

    """

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
) -> list[dict]:
    """Low level function that wraps the `ListFiles` endpoint.

    Returns a list of files from databases and rows
    based on row assigned to.

    Args:
        assigned_row_ids: ID (or human ID) or the assigned row.
        is_unassigned: Whether file is assigned to any row.

    Returns:
        A list of dictionaries, where each entry corresponds to a file. Each dictionary contains a field called `file` that corresponds conforms to a [DescribeFileResponse][src.managed_data.schema.DescribeFileResponse].

    """
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
) -> dict:
    """Low level function that wraps the `DescribeDatabaseStats` endpoint.

    Returns a dictionary of statistics about a database.


    Args:
        database_id: ID (or human ID) of the database.

    Returns:
        A dictionary that contains statistics about the database.

    """
    client = _get_default_client(client)

    return client.invoke("DescribeDatabaseStats", dict(databaseId=database_id))


@beartype
def list_mentions(
    query: str,
    *,
    client: Optional[Client] = None,
) -> dict:
    """Low level function that wraps the `ListMentions` endpoint.

    Returns a dictionary of mentions (cross references)
    of the requested object.


    Args:
        query: ID (or human ID) of row, database, workspace, or file.

    Returns:
        A dictionary that contains a field called `mentions`, which is a list of dictionaries that each refer to a row.

    """
    client = _get_default_client(client)

    return client.invoke("ListMentions", dict(query=query))


@beartype
def list_row_back_references(
    row_id: str,
    *,
    client: Optional[Client] = None,
) -> dict:
    """Low level function that wraps the `ListRowBackReferences` endpoint.

    Returns a dictionary of back references from the queried row.


    Args:
        row_id: ID (or human ID) of row.

    Returns:
        A dictionary that contains a field called `rows`, which is a list of dictionaries that each refer to a database row.

    """
    client = _get_default_client(client)

    return client.invoke("ListRowBackReferences", dict(rowId=row_id))


@beartype
def create_file_download_url(
    file_id: str,
    *,
    client: Optional[Client] = None,
) -> dict:
    """Low level function that wraps the `CreateFileDownloadUrl` endpoint.

    Returns a pre-signed URL that allows you to download a
    file.


    Args:
        file_id: ID of file.

    Returns:
        A dictionary that contains a field `downloadUrl`, that
        contains a AWS pre-signed URL.

    """

    client = _get_default_client(client)

    return client.invoke("CreateFileDownloadUrl", dict(fileId=file_id))


@beartype
def describe_file(
    file_id: str,
    *,
    client: Optional[Client] = None,
) -> dict:
    """Low level function that wraps the `DescribeFile` endpoint.

    Returns a description of file, including S3 URI, name,
    status, content length, and type.


    Args:
        file_id: ID of file.

    Returns:
        A dictionary that contains a file description, that conforms to [DescribeFileResponse][src.managed_data.schema.DescribeFileResponse].

    """

    client = _get_default_client(client)

    return client.invoke("DescribeFile", dict(fileId=file_id))


@beartype
def describe_row(
    row_id: str,
    *,
    fields: bool = False,
    client: Optional[Client] = None,
) -> dict:
    """Low level function that wraps the `DescribeRow` endpoint.

    Returns a description of a row or a database


    Args:
        row_id: ID or (human ID) or row or database.
        fields: if True, a fields item is returned in the response.

    Returns:
        A dictionary that contains a row description, that
        conforms to [DescribeRowResponse][src.managed_data.schema.DescribeRowResponse].

    """

    client = _get_default_client(client)

    data = dict(rowId=row_id, fields=fields)

    return client.invoke("DescribeRow", data)


@beartype
def list_database_rows(
    row_id: str,
    *,
    client: Optional[Client] = None,
) -> list[dict]:
    """Low level function that wraps the `ListDatabaseRows` endpoint."""

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
    """Download a file to a destination folder."""

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
    """Convert a list of human IDs to IDs or vice versa."""

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
