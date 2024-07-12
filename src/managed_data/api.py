"""The `deeporigin.managed_data.api` module contains high-level functions for
interacting with Deep Origin managed data."""

import mimetypes
import os
from pathlib import Path
from typing import Any, Optional, Union
from urllib.parse import urlparse, urlunparse

from beartype import beartype
from deeporigin import auth
from deeporigin._data_api import DeeporiginData
from deeporigin._data_api.types.rows.hierarchy_list_response import (
    Data as ListRowResponse,
)
from deeporigin.exceptions import DeepOriginException
from deeporigin.managed_data import _api
from deeporigin.managed_data.schema import (
    Cardinality,
    DataType,
    DatabaseReturnType,
    IDFormat,
    RowType,
)
from deeporigin.utils import PREFIXES


@beartype
def _get_default_client(client=None):
    """Internal function to instantiate client

    Creates and returns an authenticated client if
    not provided with one.

    Warning: Internal function
        Do not use this function

    Args:
        client: None, or a Client


    """
    if client is None:
        tokens = auth.get_tokens()
        access_token = tokens["access"]

        from deeporigin.config import get_value

        org_id = get_value()["organization_id"]

        client = DeeporiginData(
            bearer_token=access_token,
            org_id=org_id,
        )

    return client


@beartype
def list_rows(
    *,
    parent_id: Optional[str] = None,
    row_type: Optional[RowType] = None,
    parent_is_root: Optional[bool] = None,
    client=None,
) -> list[ListRowResponse]:
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

    response = client.rows.hierarchy.list(filters=filters)

    return response.data


def upload_file(
    file_path: str,
    client=None,
) -> None:
    """Upload a file to Deep Origin."""

    # attempt to guess the content type
    mime_type, _ = mimetypes.guess_type(file_path)
    content_type = mime_type if mime_type else "application/octet-stream"

    content_length = os.path.getsize(file_path)

    response = _api.create_file_upload_url(
        name=os.path.basename(file_path),
        content_type=content_type,
        content_length=content_length,
    )

    # extract pre-signed URL to upload to
    url = urlparse(response["uploadUrl"])
    url = urlunparse((url.scheme, url.netloc, url.path, "", "", ""))

    # extract params
    params = _api._parse_params_from_url(response["uploadUrl"])

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
            raise DeepOriginException(message="Error uploading file")

    return response["file"]


@beartype
def make_database_rows(
    database_id: str,
    n_rows: int = 1,
    client=None,
) -> dict:
    """Makes one or several new row(s) in a Database table

    This wraps the `EnsureRows` endpoint and sends a payload
    designed to create new row(s) in a database table.


    Args:
        database_id: ID or Human ID of the database

    Returns:
        A dictionary that conforms to a EnsureRowsResponse
    """

    if n_rows < 1:
        raise DeepOriginException(
            message=f"n_rows must be at least 1. However, n_rows was {n_rows}"
        )

    data = dict(
        databaseId=database_id,
        rows=[{"row": {}} for _ in range(n_rows)],
    )

    return _api.ensure_rows(data, client=client)


@beartype
def assign_files_to_cell(
    *,
    file_ids: list[str],
    database_id: str,
    column_id: str,
    row_id: Optional[str] = None,
    client=None,
) -> dict:
    """Assign existing file(s) to a cell

    Assign files to a cell in a database table, where the cell is identified by the database ID, row ID, and column ID. If row_id is None, a new row will be created.

    Args:
    file_id: ID of the file
    column_id: ID of the column
    row_id: ID of the row


    """

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

    return _api.ensure_rows(data, client=client)


@beartype
def upload_file_to_new_database_row(
    *,
    database_id: str,
    file_path: str,
    column_id: str,
    client=None,
):
    """Upload a file to a new row in a database.

    Upload a file to a new row in a database. This utility function
    wraps two level functions:

        - [upload_file][src.managed_data.api.upload_file]
        - [assign_files_to_cell][src.managed_data.api.assign_files_to_cell]

    Args:
        database_id: ID (or human ID) of a database.
        file_path: Path to the file to upload.
        column_id: ID (or human ID) of a column in the database.

    """
    # upload file
    response = _api.upload_file(file_path, client=client)
    file_id = response["id"]

    # assign file to column
    # we're not specifying row_id, which will create a new row
    return _api.assign_files_to_cell(
        file_ids=[file_id],
        database_id=database_id,
        column_id=column_id,
        client=client,
    )


@beartype
def get_tree(
    *,
    include_rows: bool = True,
    client=None,
) -> list[dict]:
    """Construct a tree of all workspaces, databases and rows.

    Returns a tree that contains all workspaces, databases and
    (optionally) rows. The tree is returned as a dictionary,
    and children of each object are contained in a field
    called `children`.


    Args:
        include_rows: If `True`, rows are included in the tree.

    Returns:
        A dictionary describing the tree structure of all workspaces
        and databases.

    """

    if include_rows:
        # we need to fetch everything, so use a single call
        objects = _api.list_rows(client=client)
        rows = [obj for obj in objects if obj["type"] == "row"]
        workspaces = [obj for obj in objects if obj["type"] == "workspace"]
        databases = [obj for obj in objects if obj["type"] == "database"]
    else:
        workspaces = _api.list_rows(row_type="workspace", client=client)
        databases = _api.list_rows(row_type="database", client=client)
        objects = workspaces + databases

    for obj in workspaces + databases:
        obj["children"] = []

    root_objects = [obj for obj in objects if obj["parentId"] is None]

    for root_object in root_objects:
        _add_children(root_object, workspaces)
        for workspace in workspaces:
            _add_children(workspace, workspaces + databases)

        if include_rows:
            for database in databases:
                _add_children(database, rows)

    return root_objects


@beartype
def _add_children(node: dict, objects: list[dict]) -> None:
    """Internal function to add children to a node


    Warning: Internal function
        Do not use this function.

    """
    node["children"] = [obj for obj in objects if obj["parentId"] == node["id"]]


@beartype
def get_cell_data(
    *,
    row_id: str,
    column_name: str,
    client=None,
) -> Any:
    """Extract data from a cell in a database, referenced
    by `row_id` and `column_name`.

    Returns the value in a single cell in a database.

    Warning: Caution
        This function internally calls
        [get_row_data][src.managed_data.api.get_row_data],
        so it is not efficient to write a loop to get all values
        of cells from a row. It will be faster to call
        [get_row_data][src.managed_data.api.get_row_data] directly.


    Args:
        row_id: ID (or human ID) of a row.
        column_name: Name of column.

    Returns:
        Value of that cell.

    """

    data = get_row_data(row_id, client=client)
    return data[column_name]


def set_cell_data(
    value: Any,
    *,
    database_id: str,
    row_id: str,
    column_id: str,
    client=None,
) -> Any:
    """set data in a cell to some value.

    uses the EnsureRows API endpoint"""

    # first, get the type of the column
    response = _api.describe_row(database_id, client=client)

    column = [
        col
        for col in response["cols"]
        if col["id"] == column_id or col["name"] == column_id or col["key"] == column_id
    ]

    if len(column) != 1:
        raise DeepOriginException(
            message=f"Could not find column {column_id} in database {database_id}"
        )

    column = column[0]
    if column["type"] == "select":
        if isinstance(value, list):
            pass
        elif isinstance(value, str):
            value = [value]
        elif value is None:
            value = [""]
        else:
            raise DeepOriginException(
                message="Attempting to write to a cell that is of type select. Values to be written here should be strings or lists of strings."
            )

        options = column["configSelect"]["options"]
        for item in value:
            if item not in options:
                raise DeepOriginException(
                    message=f"Expected every item to be in the options list. However, `{item}` is not in {options}"
                )

        validated_value = dict(selectedOptions=value)
    elif column["type"] == "text":
        validated_value = str(value)
    elif column["type"] == "integer":
        try:
            validated_value = int(value)
        except ValueError:
            raise DeepOriginException(
                message=f"Attempting to write to a cell that is of type integer. Value to be written here should be an integer. Instead, you attempted to write: {value}"
            )

    elif column["type"] == "float":
        try:
            validated_value = float(value)
        except ValueError:
            raise DeepOriginException(
                message=f"Attempting to write to a cell that is of type float. Value to be written here should be an float. Instead, you attempted to write: {value}"
            )

    elif column["type"] == "boolean":
        if isinstance(value, bool) or value is None:
            validated_value = value
        else:
            raise DeepOriginException(
                message=f"Attempting to write to a cell that is of type boolean. Value to be written here should be a True, False or None. Instead, you attempted to write: {value}"
            )
    else:
        raise NotImplementedError("This data type is not yet supported")

    data = {
        "databaseId": database_id,
        "rows": [
            {
                "cells": [
                    {
                        "columnId": column_id,
                        "value": validated_value,
                    }
                ],
                "row": {},
                "rowId": row_id,
            }
        ],
    }

    return _api.ensure_rows(data, client=client)


@beartype
def download(
    source: str,
    destination: str,
    *,
    include_files: bool = False,
    client=None,
) -> None:
    """Download resources from Deep Origin and save them to
    a local destination.

    Download databases, objects and other entities from
    Deep Origin managed data and save them to local disk.

    Info: Work in progress
        All features in this function have not been implemented yet.


    Args:
        source: ID (or human ID) of a resource on Deep Origin.
        destination: Path to local directory to save resources.
        include_files: if `True`, download files in database.

    """

    Path(destination).mkdir(parents=True, exist_ok=True)

    if not os.path.isdir(destination):
        raise DeepOriginException(
            message=f"{destination} should be a path to a folder."
        )

    source = source.replace(PREFIXES.DO, "")

    # first, need to determine what this is.
    if PREFIXES.FILE in source:
        # this is a file

        _api.download_file(
            file_id=source,
            destination=destination,
            client=client,
        )
        return

    # not a file, so need to determine what sort of row it is
    obj = _api.describe_row(source, client=client)
    if obj["type"] == "database":
        download_database(
            obj,
            destination,
            include_files=include_files,
            client=client,
        )
    else:
        raise NotImplementedError(
            "Downloading this type of data object has not been implemented yet"
        )


@beartype
def download_database(
    source: Union[str, dict],
    destination: str = os.getcwd(),
    *,
    include_files: bool = False,
    client=None,
) -> None:
    """Download a database and save it to a CSV file on the local disk.

    Download a database from Deep Origin managed data
    and save to local disk as a CSV file.

    Args:
        source: ID (or human ID) of a resource on Deep Origin.
        destination: Path to local directory to save resources.
        include_files: if `True`, download files in database.

    """

    if not os.path.isdir(destination):
        raise DeepOriginException(
            message=f"{destination} should be a path to a folder."
        )

    if isinstance(source, str):
        source = _api.describe_row(source, client=client)
    elif not {"hid", "id"}.issubset(set(list(source.keys()))):
        raise DeepOriginException(
            message=f"If `source` is a dictionary, expected it contain the `hid` and `id` keys. These keys were not found. Instead, the keys are: {source.keys()}"
        )

    database_id = source["id"]
    database_hid = source["hid"]
    df = get_dataframe(
        database_id,
        use_file_names=True,
        client=client,
    )

    # now download all files in the database
    if include_files:
        file_ids = df.attrs["file_ids"]

        for file_id in file_ids:
            _api.download_file(file_id, destination, client=client)

    df.to_csv(os.path.join(destination, database_hid + ".csv"))


@beartype
def get_dataframe(
    database_id: str,
    *,
    use_file_names: bool = True,
    reference_format: IDFormat = "human-id",
    return_type: DatabaseReturnType = "dataframe",
    client=None,
):
    """Generate a `pandas.DataFrame` or dictionary for a database.

    Download a database from Deep Origin managed data
    and return it as a data frame or dictionary.

    Args:
        database_id: ID (or human ID) of a database on Deep Origin.
        use_file_names: If `True`, refer to files by name rather than ID.
        reference_format: Refer to rows on Deep Origin using human IDs or system IDs.
        return_type: Whether to return a `pandas.Dataframe` or a dictionary.
    """

    # figure out the rows
    rows = _api.list_database_rows(database_id, client=client)

    # filter out template rows
    rows = [
        row for row in rows if not ("isTemplate" in row.keys() and row["isTemplate"])
    ]

    # figure out the column names and ID of the database
    response = _api.describe_row(database_id, client=client)
    assert (
        response["type"] == "database"
    ), f"Expected database_id: {database_id} to resolve to a database, but instead, it resolved to a {response['type']}"

    columns = response["cols"]
    database_id = response["id"]
    row_id = "ID"

    # make a dictionary with all data in the database
    data = dict()
    data[row_id] = []
    data["Validation Status"] = []

    # keep track of all files and references in this database
    file_ids = []
    reference_ids = []

    for column in columns:
        data[column["id"]] = []

    for row in rows:
        data[row_id].append(row["hid"])
        data["Validation Status"].append(row["validationStatus"])

        if "fields" not in row.keys():
            for column in columns:
                data[column["id"]].append(None)
            continue

        fields = row["fields"]

        for column in columns:
            value = _parse_column_value(
                column=column,
                fields=fields,
                file_ids=file_ids,
                reference_ids=reference_ids,
                use_file_names=use_file_names,
                reference_format=reference_format,
            )
            if value is None:
                data[column["id"]].append(None)
            else:
                if column["cardinality"] == "many":
                    data[column["id"]].append(value)
                else:
                    data[column["id"]].extend(value)

    if return_type == "dataframe":
        # make the dataframe

        # this import is here because we don't want to
        # import pandas unless we actually use this function
        import pandas as pd

        df = pd.DataFrame(data)
        df.attrs["file_ids"] = list(set(file_ids))
        df.attrs["reference_ids"] = list(set(reference_ids))
        df.attrs["id"] = database_id
        df.attrs["primary_key"] = row_id

        return _type_and_cleanup_dataframe(df, columns)

    else:
        # rename keys
        column_mapper = dict()
        for column in columns:
            column_mapper[column["id"]] = column["name"]
        renamed_data = data.copy()
        for key in column_mapper.keys():
            new_key = column_mapper[key]
            renamed_data[new_key] = data[key]
            renamed_data.pop(key, None)
        return renamed_data


@beartype
def _parse_column_value(
    *,
    column: dict,
    fields: list[dict],
    file_ids: list,
    reference_ids: list,
    use_file_names: bool,
    reference_format: IDFormat,
):
    """Internal function parse column values

    Warning: Internal function
        Do not use this function.
    """

    field = [field for field in fields if field["columnId"] == column["id"]]

    if len(field) == 0:
        return None

    if "value" not in field[0].keys():
        return None
    value = [field[0]["value"]]

    # special treatment for some column types
    if column["type"] == "select" and len(value) == 1:
        value = value[0]["selectedOptions"]
    elif column["type"] == "file" and len(value) == 1:
        value = value[0]["fileIds"]

        file_ids.extend(value)

        if use_file_names:
            try:
                value = [_api.describe_file(file_id)["name"] for file_id in value]
            except DeepOriginException:
                # something went wrong, ignore
                pass
    elif column["type"] == "reference" and len(value) == 1:
        value = value[0]["rowIds"]
        reference_ids.extend(value)

        if reference_format == "human-id":
            value = _api.convert_id_format(ids=value)
            value = [thing["hid"] for thing in value]

    # if there is no item
    if len(value) == 0:
        value = None

    return value


@beartype
def _type_and_cleanup_dataframe(
    df,  # pd.Dataframe, not typed to avoid pandas import
    columns: list[dict],
):
    """Internal function to type and clean a pandas dataframe

    Warning: Internal function
        Do not use this function.
    """

    # this import is here because we don't want to
    # import pandas unless we actually use this function
    import pandas as pd

    column_mapper = dict()
    for column in columns:
        column_mapper[column["id"]] = column["name"]

    for column in columns:
        col_id = column["id"]

        if column["type"] == "date":
            df[col_id] = pd.to_datetime(df[col_id])

        # special treatment for string columns
        if column["type"] in ["file", "text"]:
            df[col_id] = df[col_id].astype("string")

        if column["type"] == "boolean":
            df[col_id] = df[col_id].astype("boolean")

        # special treatment of Select columns
        if column["type"] == "select" and column["cardinality"] == "one":
            categories = column["configSelect"]["options"]
            df[col_id] = pd.Categorical(df[col_id], categories=categories)

    # rename columns
    df = df.rename(columns=column_mapper)

    # special treatment for validation status
    categories = set(df["Validation Status"])
    df["Validation Status"] = pd.Categorical(
        df["Validation Status"], categories=categories
    )

    # wipe metadata from columns
    for column in df.columns:
        df[column].attrs = dict()

    # attach metadata to columns
    for column in columns:
        df[column["name"]].attrs = column

    # add a type to Validation Status because this column
    # doesn't actually exist in the database
    df["Validation Status"].attrs = dict(type="Validation Status")

    # sort by primary key
    primary_key = df.attrs["primary_key"]
    df["numeric_id"] = df[primary_key].apply(lambda x: int(x.split("-")[1]))
    df = df.sort_values(by="numeric_id")
    df = df.drop(columns="numeric_id")
    df = df.set_index("ID")

    return df


@beartype
def get_columns(
    row_id: str,
    *,
    client=None,
) -> list[dict]:
    """Get information about the columns of a row or database.

    If `row_id` is a database, then column metadata and names
    are returned. If `row_id` is a row, then a dictionary of
    human IDs and values are returned.

    Args:
        row_id: ID (or human ID) of a row or database on Deep Origin.
    """

    response = _api.describe_row(row_id, fields=True, client=client)

    assert response["type"] in [
        "row",
        "database",
    ], "Expected row_id to resolve to a row or a database"

    if response["type"] == "database":
        # return cols info
        return response["cols"]
    else:
        return response["fields"]


@beartype
def get_row_data(
    row_id: str,
    *,
    use_column_keys: bool = False,
    client=None,
) -> dict:
    """Get the data in a row.

    Read data from a row, and return it as a dictionary, where
    the keys are column names (or keys), and the values are the values of those
    cells.

    Args:
        row_id: ID (or human ID) of a row or database on Deep Origin.
        use_column_keys: if `True`, keys of dictionary are column keys.

    Raises:
        DeepOriginException: If row_id is not a row
    """

    response = _api.describe_row(row_id, fields=True, client=client)

    if response["type"] != "row":
        raise DeepOriginException(
            message=f"Expected `row_id` to resolve to a row, instead, it resolves to a `{response['type']}`"
        )

    # ask parent for column names
    parent_response = _api.describe_row(response["parentId"], client=client)

    if parent_response["type"] != "database":
        raise DeepOriginException(
            message=f"Expected parent of `{row_id}` to resolve to a database, instead, it resolves to a `{parent_response['type']}`"
        )

    # make a dictionary from column IDs to column names
    column_name_mapper = dict()
    column_cardinality_mapper = dict()
    for col in parent_response["cols"]:
        if use_column_keys:
            column_name_mapper[col["id"]] = col["key"]
        else:
            column_name_mapper[col["id"]] = col["name"]
        column_cardinality_mapper[col["id"]] = col["cardinality"]

    # now use this to construct the required dictionary
    row_data = dict()
    for field in response["fields"]:
        column_id = field["columnId"]

        if "value" not in field:
            continue

        value = field["value"]
        if isinstance(value, dict):
            if "selectedOptions" in value.keys():
                value = value["selectedOptions"]
            elif "fileIds" in value.keys():
                value = value["fileIds"]

        if column_cardinality_mapper[column_id] == "one" and isinstance(value, list):
            value = value[0]

        row_data[column_name_mapper[column_id]] = value

    return row_data


@beartype
def merge_databases(dfs: list):
    """Merge dataframes for multiple databases into a single dataframes.

    Given a list of dataframes derived from Deep Origin databases,
    merge them into a single dataframe, resolving cross-references
    across the databases.

    Info: Work in progress
        All features in this function have not been implemented yet.


    Args:
        dfs: List of `pandas.DataFrames`.


    """

    import pandas as pd

    for df in dfs:
        assert isinstance(df, pd.DataFrame), "Expected a list of dataframes to merge"

    assert len(dfs) == 2, "For now we only support merging 2 databases"

    # make a cross reference mapper that converts
    # system IDs to column names
    cross_reference_mapper = dict()

    for df in dfs:
        cross_reference_mapper[df.attrs["id"]] = df.attrs["primary_key"]

    # rename columns that contain cross-references (foreign keys)
    # so that the pandas merge works correctly

    for df in dfs:
        column_mapper = dict()

        for column in df.columns:
            attrs = df[column].attrs

            if "referenceDatabaseRowId" in attrs.keys():
                column_mapper[column] = cross_reference_mapper[
                    attrs["referenceDatabaseRowId"]
                ]

        df.rename(columns=column_mapper, inplace=True)

    # for now we only support merging 2 DBs
    return dfs[0].merge(dfs[1], how="outer")
