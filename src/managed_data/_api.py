"""The `deeporigin.managed_data._api` module contains low-level functions for interacting
with Deep Origin's managed data API. The functions in this module
simply provide Pythonic interfaces to individual API endpoints."""

import os
from typing import Optional, Union
from urllib.parse import parse_qs, urlparse

from beartype import beartype
from deeporigin.exceptions import DeepOriginException
from deeporigin.managed_data.client import Client, DeepOriginClient
from deeporigin.managed_data.schema import (
    Cardinality,
    DataType,
    RowType,
)


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
        A dictionary that conforms to a [CreateDatabaseResponse][src.managed_data.schema.CreateDatabaseResponse].
    """

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
def convert_id_format(
    *,
    hids: Optional[Union[list[str], set[str]]] = None,
    ids: Optional[Union[list[str], set[str]]] = None,
    client: Optional[Client] = None,
) -> list[dict]:
    """Convert a list of human IDs to IDs or vice versa."""

    if hids is None and ids is None:
        raise DeepOriginException(
            message="Either `hids` or `ids` should be non-None and a list of strings"
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
