"""The `deeporigin.managed_data._api` module contains low-level functions for interacting
with Deep Origin's managed data API. The functions in this module
simply provide Pythonic interfaces to individual API endpoints."""

import mimetypes
import os
from typing import Optional, Union
from urllib.parse import parse_qs, urlparse, urlunparse

from beartype import beartype
from deeporigin.exceptions import DeepOriginException
from deeporigin.managed_data.client import Client, DeepOriginClient
from deeporigin.managed_data.schema import (
    Cardinality,
    DataType,
    RowType,
)


@beartype
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
def assign_files_to_cell(
    *,
    file_ids: list[str],
    database_id: str,
    column_id: str,
    row_id: Optional[str] = None,
    client: Optional[Client] = None,
) -> dict:
    """Assign existing file(s) to a cell

    Assign files to a cell in a database table, where the cell is identified by the database ID, row ID, and column ID. If row_id is None, a new row will be created.

    Args:
    file_id: ID of the file
    column_id: ID of the column
    row_id: ID of the row


    """
    client = _get_default_client(client)

    if row_id is None:
        data = make_database_rows(
            database_id=database_id,
            n_rows=1,
            client=client,
        )
        row_id = data["rows"][0]["id"]

    data = {
        "databaseId": database_id,
        "rows": [
            {
                "rowId": row_id,
                "row": {},
                "cells": [
                    {
                        "columnId": column_id,
                        "value": {"fileIds": file_ids},
                    },
                ],
            },
        ],
    }

    return client.invoke("EnsureRows", data)


@beartype
def make_database_rows(
    database_id: str,
    n_rows: int = 1,
    client: Optional[Client] = None,
) -> dict:
    """Makes one or several new row(s) in a Database table

    This wraps the `EnsureRows` endpoint and sends a payload
    designed to create new row(s) in a database table.


    Args:
        database_id: ID or Human ID of the database

    Returns:
        A dictionary that conforms to a [EnsureRowsResponse][src.managed_data.schema.EnsureRowsResponse]
    """

    client = _get_default_client(client)

    if n_rows < 1:
        raise DeepOriginException(
            f"n_rows must be at least 1. However, n_rows was {n_rows}"
        )

    data = dict(
        databaseId=database_id,
        rows=[{"row": {}} for _ in range(n_rows)],
    )

    return client.invoke("EnsureRows", data)


def upload_file(
    file_path: str,
    client: Optional[Client] = None,
) -> None:
    """Upload a file to Deep Origin."""

    client = _get_default_client(client)

    # attempt to guess the content type
    mime_type, _ = mimetypes.guess_type(file_path)
    content_type = mime_type if mime_type else "application/octet-stream"

    content_length = os.path.getsize(file_path)

    response = create_file_upload_url(
        name=os.path.basename(file_path),
        content_type=content_type,
        content_length=content_length,
    )

    # extract pre-signed URL to upload to
    url = urlparse(response["uploadUrl"])
    url = urlunparse((url.scheme, url.netloc, url.path, "", "", ""))

    # extract params
    params = _parse_params_from_url(response["uploadUrl"])

    headers = {
        "Accept": "application/json, text/plain, */*",
        "Connection": "keep-alive",
        "Content-Length": str(content_length),
        "Content-Type": content_type,
    }

    with open(file_path, "rb") as file:
        put_response = client.put(
            url,
            headers=headers,
            params=params,
            data=file,
        )

        if put_response.status_code != 200:
            raise DeepOriginException("Error uploading file")

    return response["file"]


@beartype
def _parse_params_from_url(url: str) -> dict:
    """utility function to extract params from a URL query

    Warning: Internal function
        Do not use this function

    Args:
        url: URL

    Returns:
        A dictionary of params
    """

    query = urlparse(url).query
    params = parse_qs(query)
    params = {key: value[0] for key, value in params.items()}
    return params


@beartype
def create_file_upload_url(
    *,
    name: str,
    content_type: str,
    content_length: int,
    client: Optional[Client] = None,
) -> dict:
    """low level function that wraps the `CreateFileUpload` endpoint.

    Creates a new file upload URL.

    Args:
        name: Name of the file
        content_type: Content type of the file
        content_length: Content length of the file

    Returns:
        A dictionary that conforms to a [CreateFileUploadResponse][src.managed_data.schema.CreateFileUploadResponse]
    """

    client = _get_default_client(client)

    data = {
        "name": name,
        "contentType": content_type,
        "contentLength": str(content_length),
    }

    return client.invoke("CreateFileUpload", data)


@beartype
def create_workspace(
    *,
    name: str,
    hid: Optional[str] = None,
    parent_id: Optional[str] = None,
    client: Optional[Client] = None,
) -> dict:
    """Low level function that wraps the `CreateWorkspace` endpoint.

    Creates a new workspace at the root level or within another workspace.

    Args:
        hid: Human ID of the workspace.
        name: Name of the workspace.
        parent_id: ID (or human ID) of the parent.

    Returns:
        A dictionary that conforms to a [CreateWorkspaceResponse][src.managed_data.schema.CreateWorkspaceResponse].
    """
    client = _get_default_client(client)

    if hid is None:
        hid = name

    data = dict(
        workspace=dict(
            hid=hid,
            name=name,
            parentId=parent_id,
        )
    )
    return client.invoke("CreateWorkspace", data)


@beartype
def create_database(
    *,
    name: str,
    hid_prefix: str,
    hid: Optional[str] = None,
    parent_id: Optional[str] = None,
    client: Optional[Client] = None,
) -> dict:
    """Low level function that wraps the `CreateDatabase` endpoint.

    Creates a new database within a workspace.

    Args:
        hid: Human ID of the database.
        name: Name of the database.
        parent_id: ID of the parent workspace.
        hid_prefix: Human ID prefix of the database. This prefix is used in every row.


    Returns:
        A dictionary that conforms to a [CreateDatabaseResponse][src.managed_data.schema.CreateDatabaseResponse]."""

    client = _get_default_client(client)

    if hid is None:
        hid = name

    data = dict(
        database=dict(
            hid=hid,
            name=name,
            parentId=parent_id,
            hidPrefix=hid_prefix,
        )
    )
    return client.invoke("CreateDatabase", data)


@beartype
def update_workspace(
    *,
    id: str,
    hid: Optional[str] = None,
    name: Optional[str] = None,
    parent_id: Optional[str] = None,
    client: Optional[Client] = None,
) -> dict:
    """Low level function that wraps the `UpdateWorkspace` endpoint.

    Updates a workspace.

    Args:
        hid: Human ID of the workspace.
        name: Name of the workspace.
        parent_id: ID (or human ID) of the parent.

    Returns:
        A dictionary that conforms to a [CreateWorkspaceResponse][src.managed_data.schema.CreateWorkspaceResponse].

    """
    client = _get_default_client(client)

    data = dict(id=id, workspace=dict())
    if hid is not None:
        data["workspace"]["hid"] = hid
    if name is not None:
        data["workspace"]["name"] = name

    return client.invoke("UpdateWorkspace", data)


@beartype
def update_database(
    *,
    id: str,
    hid_prefix: Optional[str] = None,
    hid: Optional[str] = None,
    name: Optional[str] = None,
    parent_id: Optional[str] = None,
    client: Optional[Client] = None,
) -> dict:
    """Low level function that wraps the `UpdateDatabase` endpoint.

    Updates a database.

    Args:
        hid: Human ID of the database.
        name: Name of the database.
        parent_id: ID (or human ID) of the parent.
        hid_prefix: Human ID prefix of the database. This prefix is used in every row.

    Returns:
        A dictionary that conforms to a [CreateDatabaseResponse][src.managed_data.schema.CreateDatabaseResponse].

    """
    client = _get_default_client(client)

    data = dict(id=id, database=dict())
    if hid is not None:
        data["database"]["hid"] = hid
    if name is not None:
        data["database"]["name"] = name
    if name is not None:
        data["database"]["hidPrefix"] = hid_prefix

    return client.invoke("UpdateDatabase", data)


@beartype
def delete_rows(
    row_ids: list[str],
    *,
    client: Optional[Client] = None,
) -> None:
    """Low level function that wraps the `DeleteRows` endpoint.

    Deletes rows, workspaces, or databases.

    Args:
        row_ids: A list of (system) row IDs

    Returns:
        None
    """
    client = _get_default_client(client)

    data = dict(rowIds=row_ids)
    client.invoke("DeleteRows", data)


@beartype
def add_database_column(
    *,
    database_id: str,
    name: str,
    type: DataType,
    key: Optional[str] = None,
    cardinality: Cardinality = "one",
    client: Optional[Client] = None,
) -> dict:
    """Low level function that wraps the `AddDatabaseColumn` endpoint.

    Adds a new column to a database.

    Args:
        database_id: ID of the database.
        name: Name of the column.
        type: Type of the column.
        key: Key of the column.
        cardinality: Cardinality of the column.

    Returns:
        A dictionary that conforms to a [AddDatabaseColumnResponse][src.managed_data.schema.AddDatabaseColumnResponse].
    """
    client = _get_default_client(client)

    if key is None:
        key = name

    column = dict(
        name=name,
        type=type,
        key=key,
        cardinality=cardinality,
    )

    if type == "select":
        column["configSelect"] = {
            "options": [],
            "canCreate": True,
        }

    data = dict(
        databaseId=database_id,
        column=column,
    )

    return client.invoke("AddDatabaseColumn", data)


@beartype
def delete_database_column(
    column_id: str,
    *,
    client: Optional[Client] = None,
) -> None:
    """Low level function that wraps the `DeleteDatabaseColumn` endpoint.

    Deletes a column from a database.

    Args:
        column_id: ID of the column.

    Returns:
        None
    """
    client = _get_default_client(client)

    client.invoke("DeleteDatabaseColumn", dict(columnId=column_id))


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

    client.download(url, save_path)


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
