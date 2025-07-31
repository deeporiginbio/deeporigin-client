"""The `deeporigin.data_hub.api` module contains high-level functions for
interacting with your Deep Origin data hub."""

import concurrent.futures
import os
from pathlib import Path
from typing import Any, Optional, Union
from urllib.parse import urlparse, urlunparse

from beartype import beartype
from tqdm import tqdm

from deeporigin.data_hub import _api
from deeporigin.data_hub._api import *  # noqa: F403
from deeporigin.exceptions import DeepOriginException
from deeporigin.utils import constants, network
from deeporigin.utils.core import find_last_updated_row, sha256_checksum

network.check_for_updates()


@beartype
def convert_id_format(
    *,
    hids: Optional[Union[list[str], set[str]]] = None,
    ids: Optional[Union[list[str], set[str]]] = None,
    client=None,
    _stash: bool = False,
) -> list[dict]:
    """Convert a list of human IDs to IDs or vice versa.

    Args:
        hids: List of human IDs
        ids: List of IDs (system IDs)

    """

    if hids is None and ids is None:
        raise DeepOriginException(
            message="Either the `hids` or `ids` argument should be a list of strings."
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
        _stash=_stash,
    )


@beartype
def create_workspace(
    *,
    name: str,
    hid: Optional[str] = None,
    parent_id: Optional[str] = None,
    client=None,
    _stash: bool = False,
):
    """Create a new folder (workspace) in the data hub

    A folder contains can contain other rows and databases.

    Args:
        name: Name of the folder to create
        hid: Human ID. If not specified, the name will be used
        parent_id: ID of the parent. If None, the folder is created at the root level
    """
    if hid is None:
        hid = name

    data = dict(name=name, hid=hid, parentId=parent_id)
    return _api.create_workspace(
        workspace=data,
        client=client,
        _stash=_stash,
    )


@beartype
def create_database(
    *,
    name: str,
    client=None,
    _stash: bool = False,
    parent_id: Optional[str] = None,
    hid: Optional[str] = None,
    hid_prefix: Optional[str] = None,
):
    """Create a new database in the data hub

    A database contains rows of data.

    Args:
        name: Name of the database to create
        hid: Human ID. If not specified, the name will be used
        parent_id: ID of the parent. If None, the folder is created at the root level
        hid_prefix: Human ID prefix to be used for each row. If not specified, the name will be used
    """
    if hid_prefix is None:
        hid_prefix = name

    if hid is None:
        hid = name

    data = dict(
        name=name,
        hid=hid,
        hidPrefix=hid_prefix,
        parentId=parent_id,
    )
    return _api.create_database(
        database=data,
        client=client,
        _stash=_stash,
    )


@beartype
def list_files(
    *,
    assigned_row_ids: Optional[list[str]] = None,
    is_unassigned: Optional[bool] = None,
    file_ids: Optional[list[str]] = None,
    client=None,
    _stash: bool = False,
) -> list:
    """List files, with option to filter by assigned rows, assigned status

    Returns a list of files.

    Args:
        assigned_row_ids: List of IDs of rows that files are assigned to
        is_unassigned:  If `True` only files that are unassigned will be returned

    Returns:
        A list of files, where each entry is an object that corresponds to a file on Deep Origin

    """
    filters = []

    if file_ids is not None:
        filters.append(dict(file_ids=file_ids))

    if assigned_row_ids is not None:
        filters.append(dict(assigned_row_ids=assigned_row_ids))

    if is_unassigned is not None:
        filters.append(dict(is_unassigned=is_unassigned))

    return _api.list_files(
        filters=filters,
        client=client,
        _stash=_stash,
    )


@beartype
def list_rows(
    *,
    parent_id: Optional[str] = None,
    row_type: constants.ObjectType = None,
    parent_is_root: Optional[bool] = None,
    client=None,
    _stash: bool = False,
) -> list:
    """List rows in a database or folder (workspace).

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

    return _api.list_rows(
        filters=filters,
        client=client,
        _stash=_stash,
    )


@beartype
def upload_file(
    file_path: str,
    *,
    client=None,
    _stash: bool = False,
    compute_hash: bool = True,
):
    """Upload a file to Deep Origin.

    This upload files to your Deep Origin data hub.
    To assign this file to a cell, next run [assign_files_to_cell][src.data_hub.api.assign_files_to_cell]

    Args:
        file_path: Path to the file to upload

    """

    import mimetypes

    # attempt to guess the content type
    mime_type, _ = mimetypes.guess_type(file_path)
    content_type = mime_type if mime_type else "application/octet-stream"

    content_length = os.path.getsize(file_path)

    args = dict(
        name=os.path.basename(file_path),
        content_type=content_type,
        content_length=str(content_length),
        client=client,
        _stash=_stash,
    )

    if compute_hash:
        hash = sha256_checksum(file_path)
        args["checksum_sha256"] = hash

    response = _api.create_file_upload(**args)

    # extract pre-signed URL to upload to
    url = urlparse(response.uploadUrl)
    url = urlunparse((url.scheme, url.netloc, url.path, "", "", ""))

    # extract params
    params = network._parse_params_from_url(response.uploadUrl)

    headers = {
        "Accept": "application/json, text/plain, */*",
        "Connection": "keep-alive",
        "Content-Length": str(content_length),
        "Content-Type": content_type,
    }

    if compute_hash:
        headers["x-amz-checksum-sha256"] = hash

    with open(file_path, "rb") as file:
        import requests

        put_response = requests.put(
            url,
            headers=headers,
            params=params,
            data=file,
        )

        if put_response.status_code != 200:
            raise DeepOriginException(
                title=f"❌ File could not be uploaded! [{put_response.status_code}]",
                message=f"The response was: {put_response.text}",
            )

    return response.file


@beartype
def add_database_rows(
    *,
    database_id: str,
    data: dict,
    client=None,
    _stash: bool = False,
) -> list[str]:
    """Add new data to a database.

    Use this function to add new rows, or fragments of rows, to a Deep Origin database.

    Args:
        database_id: Human ID or System ID of the database
        data: A dictionary where each key is a column name and each value is a list of values. All values should have the same length. Key names should match column names in the database.

    Returns:
        A list of row IDs

    """
    # check that dict has columns that make sense
    db = _api.describe_database(
        database_id=database_id,
        client=client,
        _stash=_stash,
    )

    col_names = [col.name for col in db.cols]

    for col in data.keys():
        if col not in col_names:
            raise DeepOriginException(
                message=f"Column `{col}` does not exist in database `{database_id}`."
            )

    # check that dict has all keys of the same length
    value_lengths = []
    for col in data.keys():
        value_lengths.append(len(data[col]))

    if len(set(value_lengths)) > 1:
        raise DeepOriginException(
            message="All rows must have the same number of values."
        )

    response = make_database_rows(
        database_id=database_id,
        n_rows=value_lengths[0],
        client=client,
        _stash=_stash,
    )

    row_ids = [row.id for row in response.rows]
    row_hids = [row.hid for row in response.rows]

    for col in data.keys():
        set_data_in_cells(
            values=data[col],
            row_ids=row_ids,
            column_id=col,
            database_id=database_id,
            columns=db.cols,
            client=client,
            _stash=_stash,
        )

    return row_hids


@beartype
def make_database_rows(
    database_id: str,
    n_rows: int = 1,
    *,
    client=None,
    _stash: bool = False,
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
            message=f"n_rows must be at least 1. However, n_rows was {n_rows}."
        )

    return _api.ensure_rows(
        client=client,
        _stash=_stash,
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
    _stash: bool = False,
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

    rows = [
        {
            "row": {},
            "cells": [
                {
                    "columnId": column_id,
                    "value": {"fileIds": file_ids},
                },
            ],
        },
    ]

    if row_id:
        rows[0]["rowId"] = row_id

    return _api.ensure_rows(
        database_id=database_id,
        rows=rows,
        client=client,
        _stash=_stash,
    )


@beartype
def upload_file_to_new_database_row(
    *,
    database_id: str,
    file_path: str,
    column_id: str,
    client=None,
    _stash: bool = False,
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
        _stash=_stash,
    )


@beartype
def get_tree(
    *,
    include_rows: bool = True,
    client=None,
    _stash: bool = False,
) -> list:
    """Construct a tree of all folders (workspaces), databases and rows.

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
        objects = list_rows(
            client=client,
            _stash=_stash,
        )
        rows = [obj for obj in objects if obj.type == "row"]
        folders = [obj for obj in objects if obj.type == "workspace"]
        databases = [obj for obj in objects if obj.type == "database"]
    else:
        folders = list_rows(
            row_type="workspace",
            client=client,
            _stash=_stash,
        )
        databases = list_rows(
            row_type="database",
            client=client,
            _stash=_stash,
        )
        objects = folders + databases

    for obj in folders + databases:
        obj["children"] = []

    # sometimes, parentId is missing (in which case it should be treated as  root object)
    root_objects = [obj for obj in objects if getattr(obj, "parentId", None) is None]

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
    node["children"] = [
        obj for obj in objects if getattr(obj, "parentId", None) == node.id
    ]


@beartype
def get_cell_data(
    *,
    row_id: str,
    column_name: str,
    client=None,
    _stash: bool = False,
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
    columns: Optional[list[dict]] = None,
    client=None,
    _stash: bool = False,
):
    """Set data in multiple cells to some value."""
    # first, get the type of the column

    if columns is None:
        # we need to get the columns and their types
        response = _api.describe_row(
            row_id=database_id,
            client=client,
            _stash=_stash,
        )

        columns = response.cols

    column = [
        col for col in columns if col["id"] == column_id or col["name"] == column_id
    ]

    if len(column) != 1:
        raise DeepOriginException(
            message=f"Column {column_id} could not be found in database {database_id}."
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
        for (row_id, validated_value) in zip(row_ids, validated_values, strict=False)
    ]

    # we cannot write more than a 1000 rows at once.
    # if there are more than 1000 rows, we need to chunk it
    max_size = 1000
    if len(rows) > max_size:
        chunks = [rows[i : i + max_size] for i in range(0, len(rows), max_size)]

        all_responses = []
        for chunk in chunks:
            all_responses.append(
                _api.ensure_rows(
                    rows=chunk,
                    client=client,
                    _stash=_stash,
                    database_id=database_id,
                )
            )

            return all_responses

    else:
        return _api.ensure_rows(
            rows=rows,
            client=client,
            _stash=_stash,
            database_id=database_id,
        )


def set_cell_data(
    value: Any,
    *,
    database_id: str,
    row_id: str,
    column_id: str,
    client=None,
    _stash: bool = False,
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
        _stash=_stash,
    )

    column = [
        col
        for col in response.cols
        if col["id"] == column_id or col["name"] == column_id
    ]

    if len(column) != 1:
        raise DeepOriginException(
            message=f"Column {column_id} could not be found in database {database_id}."
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
        _stash=_stash,
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
                message=f"Value {value} could not be written to cell {column['name']} of type select. The value of the cell must be a string or list of strings."
            )

        options = column["configSelect"]["options"]
        for item in value:
            if item not in options:
                raise DeepOriginException(
                    message=f"`{item}` is not a valid option for cell {column['name']} of type select. The valid options are {options}."
                )

        validated_value = dict(selectedOptions=value)
    elif column["type"] == "text":
        validated_value = str(value)
    elif column["type"] == "integer":
        from pandas._libs.missing import NAType

        if isinstance(value, NAType):
            validated_value = None
        else:
            try:
                validated_value = int(value)
            except ValueError as e:
                raise DeepOriginException(
                    message=f"{value} is not valid for cell {column['name']} of type integer. The value must be an integer."
                ) from e

    elif column["type"] == "float":
        from pandas._libs.missing import NAType

        if isinstance(value, NAType):
            validated_value = None
        else:
            try:
                validated_value = float(value)
            except ValueError as e:
                raise DeepOriginException(
                    message=f"{value} is not valid for cell {column['name']} of type float. The value must be a float."
                ) from e

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
                message=f"{value} is not valid for cell {column['name']} of type Boolean. The value must be a Boolean: true, false, 0, 1, yes, no, y, or n."
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
    _stash: bool = False,
) -> None:
    """Download resources from Deep Origin and save them to
    a local destination.

    Download databases, objects and other entities from
    your Deep Origin data hub and save them to local disk.

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
            message=f"Destination `{destination}` should be a path for a folder."
        )

    source = source.replace(constants.PREFIXES.DO, "")

    # first, need to determine what this is.
    if constants.PREFIXES.FILE in source:
        # this is a file

        download_files(
            file_ids=[source],
            save_to_dir=destination,
            client=client,
            _stash=_stash,
        )
        return

    # not a file, so need to determine what sort of row it is
    obj = _api.describe_row(
        row_id=source,
        client=client,
        _stash=_stash,
    )
    if obj.type == "database":
        download_database(
            obj.id,
            destination,
            include_files=include_files,
            client=client,
            _stash=_stash,
        )
    else:
        raise NotImplementedError(
            "Downloading this type of data object has not been implemented yet"
        )


@beartype
def download_database(
    source: Any,
    destination: str = os.getcwd(),
    *,
    include_files: bool = False,
    client=None,
    _stash: bool = False,
) -> None:
    """Download a database and save it to a CSV file on the local disk.

    Download a database from your Deep Origin data hub
    and save to local disk as a CSV file.

    Args:
        source: ID (or human ID) of a resource on Deep Origin.
        destination: Path to local directory to save resources.
        include_files: if `True`, download files in database.

    """

    if not os.path.isdir(destination):
        raise DeepOriginException(
            message=f"Destination `{destination}` should be a path for a folder."
        )

    if isinstance(source, str):
        source = _api.describe_row(
            row_id=source,
            client=client,
            _stash=_stash,
        )

    database_id = source.id
    database_hid = source.hid
    df = get_dataframe(
        database_id,
        use_file_names=True,
        reference_format="system-id",
        client=client,
        _stash=_stash,
    )

    # now download all files in the database
    if include_files:
        download_files(
            file_ids=df.attrs["file_ids"],
            save_to_dir=destination,
            client=client,
        )

    df.to_csv(os.path.join(destination, database_hid + ".csv"))


def _replace_with_mapper(item, mapper: dict):
    """utility function to replace items in a nested list with values generated by a mapper (a dict)."""

    if isinstance(item, list):
        return [_replace_with_mapper(sub_item, mapper) for sub_item in item]
    return mapper.get(item, item)


@beartype
def get_dataframe(
    database_id: str,
    *,
    use_file_names: bool = True,
    reference_format: constants.IDFormat = "human-id",
    return_type: constants.DatabaseReturnType = "dataframe",
    filter: Optional[dict] = None,
    client=None,
    _stash: bool = False,
):
    """Generate a `pandas.DataFrame` or dictionary for a database.

    Download a database from your Deep Origin data hub
    and return it as a data frame or dictionary.

    Args:
        database_id: ID (or human ID) of a database on Deep Origin.
        use_file_names: If `True`, refer to files by name rather than ID.
        reference_format: Refer to rows on Deep Origin using human IDs or system IDs.
        return_type: Whether to return a `pandas.Dataframe` or a dictionary.
    """

    # TODO: list_database_rows and describe_row can be called in parallel

    # figure out the column names and ID of the database
    db_row = _api.describe_row(
        row_id=database_id,
        client=client,
        _stash=_stash,
    )

    # figure out the rows
    if filter is None:
        rows = _api.list_database_rows(
            database_row_id=database_id,
            client=client,
            _stash=_stash,
        )
    else:
        # we may have to resolve column names
        column_ids = [col.id for col in db_row.cols]
        column_names = [col.name for col in db_row.cols]
        if filter.column_id not in column_ids:
            column_id = [col.id for col in db_row.cols if col.name == filter.column_id]
            if len(column_id) == 1:
                filter.column_id = column_id[0]
            else:
                raise DeepOriginException(
                    f"Filter column with ID or name: {filter.column_id} not found in database {database_id}",
                    fix=f"Valid column names are: {column_names}",
                )

        rows = _api.list_database_rows(
            database_row_id=database_id,
            client=client,
            _stash=_stash,
            filter=filter,
        )

    # filter out template rows
    rows = [
        row for row in rows if not (hasattr(row, "is_template") and row.is_template)
    ]

    if db_row.type != "database":
        raise DeepOriginException(
            f"Expected database_id: {database_id} to resolve to a database, but instead, it resolved to a {db_row.type}"
        )

    # early exit for empty DB
    if "cols" not in db_row.keys() or db_row.cols is None:
        data = dict()
        if return_type == "dataframe":
            # this import is here because we don't want to
            # import pandas unless we actually use this function
            df = _make_deeporigin_dataframe(
                data=data,
                reference_ids=None,
                file_ids=None,
                db_row=db_row,
                rows=None,
                columns=None,
            )
            return df
        else:
            return dict()

    columns = db_row.cols

    # make a dictionary with all data in the database
    data = dict()
    data["ID"] = []
    data["Validation Status"] = []

    # keep track of all files and references in this database
    reference_ids = []
    file_ids = []

    # remove notebook columns because they are not
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
        # warning: add_row_to_data mutates data, file_ids
        # and reference_ids
        add_row_to_data(
            data=data,
            row=row,
            columns=columns,
            file_ids=file_ids,
            reference_ids=reference_ids,
        )

    # make a dict that maps from file IDs to file names
    if use_file_names and len(file_ids) > 0:
        file_id_mapper = dict()

        # determine file name for every file ID in the dataframe
        files = _api.list_files(
            filters=[dict(fileIds=file_ids)],
            client=client,
            _stash=_stash,
        )
        for file in files:
            file_id_mapper[file.file.id] = file.file.name

        for column in columns:
            if column["type"] == "file":
                inputs = data[column["id"]]

                data[column["id"]] = [
                    _replace_with_mapper(item, file_id_mapper) for item in inputs
                ]

    # make a dict that maps row system IDs to human IDs
    if reference_format == "human-id":
        ref_id_mapper = dict()
        conversions = convert_id_format(ids=reference_ids)
        for conversion in conversions:
            ref_id_mapper[conversion.id] = conversion.hid

        for column in columns:
            if column["type"] == "reference":
                inputs = data[column["id"]]

                data[column["id"]] = [
                    _replace_with_mapper(item, ref_id_mapper) for item in inputs
                ]

    if return_type == "dataframe":
        # make the dataframe

        df = _make_deeporigin_dataframe(
            data=data,
            reference_ids=reference_ids,
            file_ids=file_ids,
            db_row=db_row,
            rows=rows,
            columns=columns,
        )
        return df

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


def _make_deeporigin_dataframe(
    *,
    data: dict,
    reference_ids: Optional[list],
    file_ids: Optional[list],
    db_row: dict,
    columns: Optional[list],
    rows: Optional[list],
):
    # this import is here because we don't want to
    # import pandas unless we actually use this function
    from deeporigin.data_hub.dataframe import DataFrame

    df = DataFrame(data)
    if reference_ids is not None:
        df.attrs["reference_ids"] = list(set(reference_ids))
        df.attrs["reference_ids"].sort()

    if file_ids is not None:
        df.attrs["file_ids"] = list(set(file_ids))
        df.attrs["file_ids"].sort()

    df.attrs["id"] = db_row.id
    df.attrs["metadata"] = dict(db_row)

    if columns is not None:
        df = _type_and_cleanup_dataframe(df, columns)

    # find last updated row for pretty printing
    if len(df) > 0:
        df.attrs["last_updated_row"] = find_last_updated_row(rows)
    else:
        df.attrs["last_updated_row"] = db_row

    df._deep_origin_out_of_sync = False
    df._modified_columns = dict()
    return df


@beartype
def download_files(
    *,
    files: Optional[list[dict]] = None,
    file_ids: Optional[list[str]] = None,
    save_to_dir: Path | str = Path("."),
    use_file_names: bool = True,
    client=None,
    _stash: bool = False,
) -> None:
    """download multiple files in parallel to local disk

    Args:
        files: list of files to download. These can be a list of file_ids or a list of files as returned by api.list_files
        save_to_dir: directory to save files to on local computer
        use_file_names: If `True`, refer to files by name rather than ID.
    """

    if not os.path.isdir(save_to_dir):
        raise DeepOriginException(
            message=f"Destination `{save_to_dir}` should be a path for a folder."
        )

    if files is None and file_ids is None:
        # nothing provided, download everything
        files = list_files(client=client, _stash=_stash)
    elif files is not None and file_ids is None:
        # list of files provided
        pass
    elif files is None and file_ids is not None:
        # list of file IDs provided
        files = list_files(
            file_ids=file_ids,
            client=client,
            _stash=_stash,
        )

    else:
        raise DeepOriginException("Only one of `files` or `file_ids` can be provided")

    if isinstance(save_to_dir, str):
        save_to_dir = Path(save_to_dir)

    file_ids = [item.file.id for item in files]

    if use_file_names:
        save_paths = [save_to_dir / item.file.name for item in files]
    else:
        save_paths = [
            save_to_dir / item.file.id.replace("_file:", "") for item in files
        ]

    with concurrent.futures.ThreadPoolExecutor() as executor:
        futures = []
        for file_id, save_path in zip(file_ids, save_paths, strict=False):
            futures.append(
                executor.submit(
                    lambda file_id, save_path: network.download_sync(
                        _api.create_file_download_url(
                            file_id=file_id, client=client
                        ).downloadUrl,
                        save_path,
                    ),
                    file_id,
                    save_path,
                )
            )

        for _ in tqdm(
            concurrent.futures.as_completed(futures),
            total=len(futures),
            desc="Downloading files from Deep Origin",
        ):
            pass


@beartype
def add_row_to_data(
    *,
    data: dict,
    row: dict,
    columns: list,
    file_ids: list,
    reference_ids: list,
):
    """utility function to combine data from a row into a dataframe"""
    row_data = row_to_dict(
        row,
        file_ids=file_ids,
        reference_ids=reference_ids,
    )
    if row_data is None:
        for column in columns:
            col_id = column["id"]
            data[col_id].append(None)

        return

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


@beartype
def row_to_dict(
    row: dict,
    *,
    file_ids: Optional[list] = None,
    reference_ids: Optional[list] = None,
) -> dict:
    """convert a database row (as returned by api.list_database_rows) to a dictionary where keys are column IDs and values are the values in the row

    Danger: This function mutates inputs
        This function mutates file_ids and reference_ids

    Args:
        row: database row (as returned by api.list_database_rows)
        file_ids: list of file IDs, will be mutated in-place
        reference_ids: list of reference IDs, will be mutated in-place

    Returns:
        dict
    """

    if file_ids is None:
        file_ids = []
    if reference_ids is None:
        reference_ids = []

    values = {"ID": row.hid, "Validation Status": row.validationStatus}

    if "fields" not in row.keys() or row.fields is None:
        return values

    fields = row.fields

    for field in fields:
        if "systemType" in field.keys() and field.systemType == "bodyDocument":
            continue
        if getattr(field, "value", None) is None:
            value = None
        elif field.type in ["float", "integer", "boolean"]:
            value = field.value
        elif field.type == "select":
            value = field.value.selectedOptions

        elif field.type == "reference":
            value = field.value.rowIds
            reference_ids.extend(value)
        elif field.type == "file":
            value = field.value.fileIds
            file_ids.extend(value)
        elif field.type == "expression":
            try:
                value = field.value.expression
            except Exception:
                value = None

        elif field.type == "user":
            user_ids = field.value.userDrns

            value = [user_id for user_id in user_ids]
        else:
            value = field.value
        values[field.columnId] = value
    return values


@beartype
def _type_and_cleanup_dataframe(
    df,  # Dataframe, not typed to avoid pandas import
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
        elif column["type"] in ["file", "text"]:
            df[col_id] = df[col_id].astype("string")

        elif column["type"] == "boolean":
            df[col_id] = df[col_id].astype("boolean")

        # special treatment of Select columns
        elif column["type"] == "select" and column["cardinality"] == "one":
            categories = column["configSelect"]["options"]
            df[col_id] = pd.Categorical(df[col_id], categories=categories)

        elif column["type"] == "integer":
            df[col_id] = df[col_id].astype("Int64")

        elif column["type"] == "float":
            df[col_id] = df[col_id].astype("Float64")

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
    _stash: bool = False,
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
        _stash=_stash,
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
def get_notebook(
    row_id: str,
    *,
    client=None,
    _stash: bool = False,
) -> list:
    """Get the notebook of a row, if it exists

    Args:
        row_id: ID (or human ID) of a row on Deep Origin.

    Returns:
        The notebook of the row, returned as a list
        of blocks

    """

    response = _api.describe_row(
        row_id=row_id,
        fields=True,
        client=client,
        _stash=_stash,
    )

    for field in response.fields:
        if field["systemType"] == "bodyDocument":
            return field["value"]["topLevelBlocks"]


@beartype
def get_row_data(
    row_id: str,
    *,
    use_column_keys: bool = False,
    client=None,
    _stash: bool = False,
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
        _stash=_stash,
    )

    if response.type != "row":
        raise DeepOriginException(
            message=f"The row could not be retrieved because `{row_id}` is a `{response.type}`. row_id must be an ID for a row."
        )

    # ask parent for column names
    parent_response = _api.describe_row(
        row_id=response.parentId,
        client=client,
        _stash=_stash,
    )

    if parent_response.type != "database":
        raise DeepOriginException(
            message=f"The row could not be retrieved because the parent of `{row_id}` is a `{parent_response.type}`. The parent of row_id must be a database."
        )

    # make a dictionary from column IDs to column names
    column_name_mapper = dict()
    column_cardinality_mapper = dict()
    for col in parent_response.cols:
        column_name_mapper[col["id"]] = col["name"]
        column_cardinality_mapper[col["id"]] = col["cardinality"]

    # now use this to construct the required dictionary

    row_data = dict()
    for col in parent_response.cols:
        if "systemType" in col.keys() and col["systemType"] == "bodyDocument":
            continue
        row_data[col["name"]] = None
    if not hasattr(response, "fields"):
        return row_data
    for field in response.fields:
        if "systemType" in field.keys() and field["systemType"] == "bodyDocument":
            continue

        column_id = field["columnId"]

        if "value" not in field:
            value = None
        else:
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
def add_database_column(
    *,
    database_id: str,
    type: constants.DataType,
    name: str,
    cardinality: constants.Cardinality = "one",
    required: bool = False,
    client=None,
    _stash: bool = False,
):
    """Add a column to a database.

    Args:
        database_id: ID (or human ID) of a database on Deep Origin.
        type: type of the column. Should be one of [DataType](types.md#src.utils.constants.DataType)
        name: name of the column
        cardinality: cardinality of the column. Specifies whether cells in this column can contain or many items. Should be one of "one" or "many"
        required: whether the column is required. If True, cells in this column cannot be empty


    """
    column = dict(
        name=name,
        type=type,
        isRequired=required,
        cardinality=cardinality,
    )

    response = _api.add_database_column(
        column=column,
        database_id=database_id,
        client=client,
        _stash=_stash,
    )

    return response


@beartype
def add_smiles_column(
    *,
    database_id: str,
    name: str,
    client=None,
    _stash: bool = False,
):
    """Add a SMILES column, with a inline 2D viewer configured

    Args:
        database_id: ID (or human ID) of a database on Deep Origin.
        name: name of the column


    """
    column = {
        "cardinality": "one",
        "cellJsonSchema": {"format": "smiles"},
        "inlineViewer": "molecule2d",
        "name": name,
        "type": "text",
    }

    return _api.add_database_column(
        column=column,
        database_id=database_id,
    )
