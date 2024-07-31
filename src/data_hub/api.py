"""The `deeporigin.data_hub.api` module contains high-level functions for
interacting with the Deep Origin data hub."""

import mimetypes
import os
from pathlib import Path
from typing import Any, Optional, Union
from urllib.parse import urlparse, urlunparse

from beartype import beartype

# this import is to allow us to use functions
# not marked in __all__ in _api
from deeporigin.data_hub import _api

# this import is to make sure that all simply-wrapped
# functions in _api are available in api (this module)
from deeporigin.data_hub._api import *  # noqa: F403
from deeporigin.exceptions import DeepOriginException
from deeporigin.utils import (
    PREFIXES,
    DatabaseReturnType,
    IDFormat,
    RowType,
    _parse_params_from_url,
    download_sync,
)


def ensure_client(func):
    """decorator to make sure that the client is configured"""

    def wrapper(*args, **kwargs):
        if "client" not in kwargs or kwargs["client"] is None:
            kwargs["client"] = _api._get_default_client()
        return func(*args, **kwargs)

    return wrapper


@beartype
@ensure_client
def convert_id_format(
    *,
    hids: Optional[Union[list[str], set[str]]] = None,
    ids: Optional[Union[list[str], set[str]]] = None,
    client=None,
) -> list[dict]:
    """Convert a list of human IDs to IDs or vice versa.

    Args:
        hids: List of human IDs
        ids: List of IDs (system IDs)

    """

    if hids is None and ids is None:
        raise DeepOriginException(
            message="Either `hids` or `ids` should be non-None and a list of strings"
        )

    conversions = []

    if hids is not None:
        for hid in hids:
            conversions.append(dict(hid=hid))

    if ids is not None:
        for sid in ids:
            conversions.append(dict(id=sid))

    return _api.convert_id_format(
        conversions=conversions,
        client=client,
    )


@beartype
@ensure_client
def list_rows(
    *,
    parent_id: Optional[str] = None,
    row_type: RowType = None,
    parent_is_root: Optional[bool] = None,
    client=None,
) -> list:
    """List rows in a database or folder.

    Returns a list of rows from folders and databases,
    based on the parent, row type, or whether the parent is the root.

    Args:
        parent_id: ID (or human ID) or the parent.
        row_type: One of `row`, `folder`, or `database`.
        parent_is_root: If `True` only rows that are children of the root will be returned.

    Returns:
        A list of objects, where each entry corresponds to a row.

    """
    filters = []

    if parent_is_root is not None:
        filters.append(dict(parent=dict(is_root=parent_is_root)))

    if parent_id:
        filters.append(dict(parent=dict(id=parent_id)))

    if row_type:
        filters.append(dict(row_type=row_type))

    return _api.list_rows(filters=filters, client=client)


@beartype
@ensure_client
def download_file(
    file_id: str,
    *,
    destination: str = os.getcwd(),
    client=None,
) -> None:
    """Download a file to a destination folder.

    Download a file synchronously from Deep Origin
    to folder on the local file system.

    Args:
        file_id: ID of the file on Deep Origin
        destination: Path to the destination folder

    """

    if not os.path.isdir(destination):
        raise DeepOriginException(
            message=f"{destination} should be a path to a folder."
        )

    file_name = _api.describe_file(file_id=file_id, client=client).name

    url = _api.create_file_download_url(file_id=file_id, client=client).download_url

    save_path = os.path.join(destination, file_name)

    download_sync(url, save_path)


@beartype
@ensure_client
def upload_file(
    file_path: str,
    client=None,
) -> None:
    """Upload a file to Deep Origin.

    This uploads to the "staging area" of the Deep Origin data hub.
    To assign this file to a cell, use [assign_files_to_cell][src.data_hub.api.assign_files_to_cell]

    Args:
        file_path: Path to the file to upload

    """

    # attempt to guess the content type
    mime_type, _ = mimetypes.guess_type(file_path)
    content_type = mime_type if mime_type else "application/octet-stream"

    content_length = os.path.getsize(file_path)

    response = _api.create_file_upload(
        name=os.path.basename(file_path),
        content_type=content_type,
        content_length=str(content_length),
        client=client,
    )

    # extract pre-signed URL to upload to
    url = urlparse(response.upload_url)
    url = urlunparse((url.scheme, url.netloc, url.path, "", "", ""))

    # extract params
    params = _parse_params_from_url(response.upload_url)

    headers = {
        "Accept": "application/json, text/plain, */*",
        "Connection": "keep-alive",
        "Content-Length": str(content_length),
        "Content-Type": content_type,
    }

    with open(file_path, "rb") as file:
        import requests

        put_response = requests.put(
            url,
            headers=headers,
            params=params,
            data=file,
        )

        if put_response.status_code != 200:
            raise DeepOriginException(message="Error uploading file")

    return response.file


@beartype
@ensure_client
def make_database_rows(
    database_id: str,
    n_rows: int = 1,
    client=None,
) -> dict:
    """Makes one or several new row(s) in a database table



    Args:
        database_id: ID or Human ID of the database
        n_rows: Number of rows to create. Must be an integer greater than 0

    Returns:
        A dictionary that conforms to a EnsureRowsResponse
    """

    if n_rows < 1:
        raise DeepOriginException(
            message=f"n_rows must be at least 1. However, n_rows was {n_rows}"
        )

    return _api.ensure_rows(
        client=client,
        rows=[{"row": {}} for _ in range(n_rows)],
        database_id=database_id,
    )


@beartype
def assign_files_to_cell(
    *,
    file_ids: list[str],
    database_id: str,
    column_id: str,
    row_id: Optional[str] = None,
    client=None,
):
    """Assign existing file(s) to a cell

    Assign files to a cell in a database table, where the cell is identified by the database ID, row ID, and column ID.

    If row_id is `None`, a new row will be created.

    Args:
        file_ids: ID of the file
        database_id: ID of database to assign to
        column_id: ID of the column
        row_id: ID of the row


    """

    if row_id is None:
        data = make_database_rows(
            database_id=database_id,
            n_rows=1,
            client=client,
        )
        row_id = data.rows[0].id

    rows = [
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
    ]

    return _api.ensure_rows(
        database_id=database_id,
        rows=rows,
        client=client,
    )


@beartype
@ensure_client
def upload_file_to_new_database_row(
    *,
    database_id: str,
    file_path: str,
    column_id: str,
    client=None,
):
    """Upload a file to a new row in a database.

    Upload a file to a new row in a database. This utility
    function wraps two other functions:

    - [upload_file][src.data_hub.api.upload_file]
    - [assign_files_to_cell][src.data_hub.api.assign_files_to_cell]

    Args:
        database_id: ID (or human ID) of a database.
        file_path: Path to the file to upload.
        column_id: ID (or human ID) of a column in the database.

    """
    # upload file
    response = upload_file(file_path, client=client)
    file_id = response.id

    # assign file to column
    # we're not specifying row_id, which will create a new row
    return assign_files_to_cell(
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
) -> list:
    """Construct a tree of all folders, databases and rows.

    Returns a tree that contains all folders, databases and
    (optionally) rows. The tree is returned as a dictionary,
    and children of each object are contained in a field
    called `children`.


    Args:
        include_rows: If `True`, rows are included in the tree.

    Returns:
        A dictionary describing the tree structure of all folders
        and databases.

    """

    if include_rows:
        # we need to fetch everything, so use a single call
        objects = list_rows(client=client)
        rows = [obj.dict() for obj in objects if obj.type == "row"]
        folders = [obj for obj in objects if obj.type == "workspace"]
        databases = [obj for obj in objects if obj.type == "database"]
    else:
        folders = list_rows(row_type="workspace", client=client)
        databases = list_rows(row_type="database", client=client)
        objects = folders + databases

    # convert everything into a dict
    folders = [obj.dict() for obj in folders]
    databases = [obj.dict() for obj in databases]

    for obj in folders + databases:
        obj["children"] = []

    root_objects = [obj.dict() for obj in objects if obj.parent_id is None]

    for root_object in root_objects:
        _add_children(root_object, folders)
        for folder in folders:
            _add_children(folder, folders + databases)

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
    node["children"] = [obj for obj in objects if obj["parent_id"] == node["id"]]


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
        [get_row_data][src.data_hub.api.get_row_data],
        so it is not efficient to write a loop to get all values
        of cells from a row. It will be faster to call
        [get_row_data][src.data_hub.api.get_row_data] directly.


    Args:
        row_id: ID (or human ID) of a row.
        column_name: Name of column.

    Returns:
        Value of that cell.

    """

    data = get_row_data(row_id, client=client)
    return data[column_name]


def set_data_in_cells(
    *,
    values: list,
    row_ids: list[str],
    column_id: str,
    database_id: str,
    client=None,
):
    """Set data in multiple cells to some value."""
    # first, get the type of the column
    response = _api.describe_row(
        row_id=database_id,
        client=client,
    )

    column = [
        col
        for col in response.cols
        if col["id"] == column_id or col["name"] == column_id or col["key"] == column_id
    ]

    if len(column) != 1:
        raise DeepOriginException(
            message=f"Could not find column {column_id} in database {database_id}"
        )

    column = column[0]
    validated_values = [
        _validate_value_for_column(value=value, column=column) for value in values
    ]

    rows = [
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
        for (row_id, validated_value) in zip(row_ids, validated_values)
    ]

    return _api.ensure_rows(
        rows=rows,
        client=client,
        database_id=database_id,
    )


def set_cell_data(
    value: Any,
    *,
    database_id: str,
    row_id: str,
    column_id: str,
    client=None,
) -> Any:
    """Set data in a cell to some value.


    Args:
        value: Value to set in the cell
        database_id: ID (or human ID) of a database
        row_id: ID (or human ID) of a row
        column_id: ID (or human ID) of a column


    """

    # first, get the type of the column
    response = _api.describe_row(
        row_id=database_id,
        client=client,
    )

    column = [
        col
        for col in response.cols
        if col["id"] == column_id or col["name"] == column_id or col["key"] == column_id
    ]

    if len(column) != 1:
        raise DeepOriginException(
            message=f"Could not find column {column_id} in database {database_id}"
        )

    column = column[0]
    validated_value = _validate_value_for_column(
        value=value,
        column=column,
    )

    rows = [
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
    ]

    return _api.ensure_rows(
        rows=rows,
        client=client,
        database_id=database_id,
    )


def _validate_value_for_column(*, column: dict, value: Any):
    """helper function that validates a value for a column,
    so that it can be written to the database."""

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
        elif isinstance(value, str) and value in [
            "True",
            "true",
            "1",
            "Yes",
            "yes",
            "Y",
            "y",
        ]:
            validated_value = True
        elif isinstance(value, str) and value in [
            "False",
            "false",
            "0",
            "No",
            "no",
            "N",
            "n",
        ]:
            validated_value = False
        else:
            raise DeepOriginException(
                message=f"Attempting to write to a cell that is of type Boolean. Value to be written here should be a True, False or None. Instead, you attempted to write: {value}"
            )
    else:
        raise NotImplementedError("This data type is not yet supported")

    return validated_value


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
    the Deep Origin data hub and save them to local disk.

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

        download_file(
            file_id=source,
            destination=destination,
            client=client,
        )
        return

    # not a file, so need to determine what sort of row it is
    obj = _api.describe_row(
        row_id=source,
        client=client,
    )
    if obj.type == "database":
        download_database(
            obj.id,
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
    source: str,
    destination: str = os.getcwd(),
    *,
    include_files: bool = False,
    client=None,
) -> None:
    """Download a database and save it to a CSV file on the local disk.

    Download a database from the Deep Origin data hub
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
        source = _api.describe_row(
            row_id=source,
            client=client,
        )

    database_id = source.id
    database_hid = source.hid
    df = get_dataframe(
        database_id,
        use_file_names=True,
        client=client,
    )

    # now download all files in the database
    if include_files:
        file_ids = df.attrs["file_ids"]

        for file_id in file_ids:
            download_file(file_id, destination, client=client)

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

    Download a database from the Deep Origin data hub
    and return it as a data frame or dictionary.

    Args:
        database_id: ID (or human ID) of a database on Deep Origin.
        use_file_names: If `True`, refer to files by name rather than ID.
        reference_format: Refer to rows on Deep Origin using human IDs or system IDs.
        return_type: Whether to return a `pandas.Dataframe` or a dictionary.
    """

    # figure out the rows
    rows = _api.list_database_rows(
        database_row_id=database_id,
        client=client,
    )

    # filter out template rows
    rows = [
        row for row in rows if not (hasattr(row, "is_template") and row.is_template)
    ]

    # figure out the column names and ID of the database
    db_row = _api.describe_row(
        row_id=database_id,
        client=client,
    )
    assert (
        db_row.type == "database"
    ), f"Expected database_id: {database_id} to resolve to a database, but instead, it resolved to a {db_row.type}"

    columns = db_row.cols
    database_id = db_row.id

    # make a dictionary with all data in the database
    data = dict()
    data["ID"] = []
    data["Validation Status"] = []

    # keep track of all files and references in this database
    file_ids = []
    reference_ids = []

    # remove body document columns because they are not
    # shown in the UI as columns
    columns = [
        col
        for col in columns
        if "systemType" not in col.keys() or col["systemType"] != "bodyDocument"
    ]

    # create empty lists for each column
    for column in columns:
        data[column["id"]] = []

    for row in rows:
        add_row_to_data(
            data=data,
            row=row,
            columns=columns,
            use_file_names=use_file_names,
        )

    if return_type == "dataframe":
        # make the dataframe

        # this import is here because we don't want to
        # import pandas unless we actually use this function
        import pandas as pd

        df = pd.DataFrame(data)
        df.attrs["file_ids"] = list(set(file_ids))
        df.attrs["reference_ids"] = list(set(reference_ids))
        df.attrs["id"] = database_id

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
    fields: Optional[list],
    file_ids: list,
    reference_ids: list,
    use_file_names: bool,
    reference_format: IDFormat,
):
    """Internal function parse column values

    Warning: Internal function
        Do not use this function.
    """

    if fields is None:
        return None

    field = [field for field in fields if field.column_id == column["id"]]

    if len(field) == 0:
        return None

    if not hasattr(field[0], "value"):
        return None
    value = [field[0].value]

    # special treatment for some column types
    if column["type"] == "select" and len(value) == 1 and value[0] is not None:
        value = value[0].selected_options
    elif column["type"] == "file" and len(value) == 1 and value[0] is not None:
        value = value[0].file_ids

        file_ids.extend(value)

        if use_file_names:
            try:
                value = [_api.describe_file(file_id=file_id).name for file_id in value]
            except DeepOriginException:
                # something went wrong, ignore
                pass
    elif column["type"] == "reference" and len(value) == 1 and value[0] is not None:
        value = value[0].row_ids
        reference_ids.extend(value)

        if reference_format == "human-id":
            value = convert_id_format(ids=value)
            value = [thing.hid for thing in value]

    if len(value) == 0:
        value = None

    return value


def add_row_to_data(*, data: dict, row, columns: list, use_file_names: bool = True):
    """utility function to combine data from a row into a dataframe"""
    row_data = _row_to_dict(row, use_file_names=use_file_names)
    data["ID"].append(row_data["ID"])
    data["Validation Status"].append(row_data["Validation Status"])

    for column in columns:
        col_id = column["id"]
        if col_id in row_data.keys():
            value = row_data[col_id]
            if column["cardinality"] == "one" and isinstance(value, list):
                value = value[0]
            data[col_id].append(value)
        else:
            data[col_id].append(None)


def _row_to_dict(row, *, use_file_names: bool = True):
    """utility function to convert a row to a dictionary"""
    fields = row.fields

    values = {"ID": row.hid, "Validation Status": row.validation_status}
    if fields is None:
        return values
    for field in fields:
        if field.system_type == "bodyDocument":
            continue
        if field.value is None:
            value = None
        elif field.type in ["float", "int", "boolean"]:
            value = field.value
        elif field.type == "select":
            value = field.value.selected_options

        elif field.type == "reference":
            value = field.value.row_ids

        elif field.type == "file":
            value = field.value.file_ids
            if use_file_names and value is not None:
                try:
                    value = [
                        _api.describe_file(file_id=file_id).name for file_id in value
                    ]
                except DeepOriginException:
                    # something went wrong, ignore
                    pass
        else:
            value = field.value
        values[field.column_id] = value
    return values


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
    df["numeric_id"] = df["ID"].apply(lambda x: int(x.split("-")[-1]))
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

    response = _api.describe_row(
        row_id=row_id,
        fields=True,
        client=client,
    )

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

    response = _api.describe_row(
        row_id=row_id,
        fields=True,
        client=client,
    )

    if response.type != "row":
        raise DeepOriginException(
            message=f"Expected `row_id` to resolve to a row, instead, it resolves to a `{response.type}`"
        )

    # ask parent for column names
    parent_response = _api.describe_row(
        row_id=response.parent_id,
        client=client,
    )

    if parent_response.type != "database":
        raise DeepOriginException(
            message=f"Expected parent of `{row_id}` to resolve to a database, instead, it resolves to a `{parent_response.type}`"
        )

    # make a dictionary from column IDs to column names
    column_name_mapper = dict()
    column_cardinality_mapper = dict()
    for col in parent_response.cols:
        if use_column_keys:
            column_name_mapper[col["id"]] = col["key"]
        else:
            column_name_mapper[col["id"]] = col["name"]
        column_cardinality_mapper[col["id"]] = col["cardinality"]

    # now use this to construct the required dictionary
    row_data = dict()
    if not hasattr(response, "fields"):
        return row_data
    for field in response.fields:
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
