"""this module contains high-level functions used to
interact with the data product"""

import os
from typing import Optional, Union

import pandas as pd
from beartype import beartype
from deeporigin.exceptions import DeepOriginException
from deeporigin.managed_data._api import (
    _invoke,
    describe_file,
    describe_row,
    download_file,
    list_database_rows,
    list_rows,
)
from deeporigin.utils import PREFIX


@beartype
def get_cell_data(*, row_id: str, column_name: str):
    """extract data from a cell in a database, referenced
    by row_id and column_name"""

    data = get_row_data(row_id)
    return data[column_name]


@beartype
def download(
    source: str,
    destination: str,
    *,
    include_files: bool = False,
):
    """download resources from source and save to destination"""

    if not os.path.isdir(destination):
        raise DeepOriginException(f"{destination} should be a path to a folder.")

    source = source.replace(PREFIX, "")

    # first, need to determine what this is.
    obj = describe_row(source)
    if obj["type"] == "database":
        download_database(
            obj,
            destination,
            include_files=include_files,
        )
    else:
        raise NotImplementedError(
            "Downloading this type of data object has not been implemented yet"
        )


@beartype
def download_database(
    source: Union[str, dict],
    destination: str,
    *,
    include_files: bool = False,
) -> None:
    """download a database and save it to a CSV

    if include_files is True, files are also downloaded
        and saved to the same folder as the CSV is in
    """

    if not os.path.isdir(destination):
        raise DeepOriginException(f"{destination} should be a path to a folder.")

    if isinstance(source, str):
        source = describe_row(source)

    database_id = source["id"]
    database_hid = source["hid"]
    df = get_dataframe(database_id, use_file_names=True)

    # now download all files in the database
    if include_files:
        file_ids = df.attrs["file_ids"]

        for file_id in file_ids:
            download_file(file_id, destination)

    df.to_csv(os.path.join(destination, database_hid + ".csv"))


@beartype
def upload(source: str, destination: str) -> None:
    """upload resources from source to destination"""

    pass


@beartype
def get_children(
    objects: Optional[Union[list[dict], str]] = None,
) -> list[dict]:
    """recursively find all workspaces, databases, rows"""
    if objects is None:
        objects = get_workspaces()
    elif isinstance(objects, str):
        # need to convert a string to a object
        obj = describe_row(objects)
        obj = {key: obj[key] for key in ["id", "name", "type", "hid"]}
        objects = [obj]

    for obj in objects:
        children = list_rows(obj["id"])
        if len(children) == 0:
            continue
        obj["children"] = get_children(children)

    return objects


@beartype
def get_workspaces() -> list[dict]:
    """list workspaces in root"""

    data = dict(filters=[dict(parent=dict(isRoot=True))])
    objects = _invoke("ListRows", data)
    return [obj for obj in objects if obj["type"] == "workspace"]


@beartype
def get_dataframe(
    database_id: str,
    *,
    use_file_names: bool = True,
) -> pd.DataFrame:
    """return a dataframe of all rows in a database

    Arguments:
    -database_id: ID of database
    -use_file_name: if True, cells containing files will
        contain the name of the original uploaded file.
        if False, fileIDs are used

    """

    # figure out the rows
    rows = list_database_rows(database_id)

    # figure out the columns
    columns = get_columns(database_id)

    # make a dictionary with all data in the database
    data = dict()
    data["validation_status"] = []
    data["row"] = []

    # keep track of all files in this database
    file_ids = []

    for column in columns:
        data[column["id"]] = []

    for row in rows:
        data["row"].append(row["hid"])
        data["validation_status"].append(row["validationStatus"])

        fields = row["fields"]

        for column in columns:
            value = _parse_column_value(
                column,
                fields,
                file_ids,
                use_file_names,
            )
            if column["cardinality"] == "many":
                data[column["id"]].append(value)
            else:
                data[column["id"]].extend(value)

    # make the dataframe
    df = pd.DataFrame(data)
    df.attrs["file_ids"] = file_ids

    return _type_and_cleanup_dataframe(df, columns)


@beartype
def _parse_column_value(
    column: dict,
    fields: list[dict],
    file_ids: list,
    use_file_names: bool,
):
    """utility function to parse value of a column"""

    value = [field["value"] for field in fields if field["columnId"] == column["id"]]

    # special treatment for some column types
    if column["type"] == "select" and len(value) == 1:
        value = value[0]["selectedOptions"]
    elif column["type"] == "file" and len(value) == 1:
        value = value[0]["fileIds"]

        file_ids.extend(value)

        if use_file_names:
            value = [describe_file(file_id)["name"] for file_id in value]

    # if there is no item
    if len(value) == 0:
        value = [None]

    return value


@beartype
def _type_and_cleanup_dataframe(
    df: pd.DataFrame,
    columns: list[dict],
) -> pd.DataFrame:
    """utility function to clean up the dataframe and
    make it more usable"""

    column_mapper = dict()
    for column in columns:
        column_mapper[column["id"]] = column["name"]

    # special treatment of Select columns
    for column in columns:
        col_id = column["id"]

        if column["type"] == "select" and column["cardinality"] == "one":
            categories = column["configSelect"]["options"]
            df[col_id] = pd.Categorical(df[col_id], categories=categories)

    # rename columns
    df = df.rename(columns=column_mapper)
    df = df.set_index("row")

    # attach metadata to columns
    for column in columns:
        df[column["name"]].attrs = column

    # add a type to validation_status because this column
    # doesn't actually exist in the database
    df["validation_status"].attrs = dict(type="validation_status")
    return df


@beartype
def get_columns(row_id: str) -> list[dict]:
    """return column information.

    if row_id is a database, then column metadata and names are returned.
    if row_id is a row, then a dictionary of hids and values are returned"""

    response = describe_row(row_id, fields=True)

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
def get_row_data(row_id: str) -> dict:
    """name needs improving? this returns fields in this
    row as a dictionary, where keys are HIDs

    if row_id is not a row, an error is raised.
    """

    response = describe_row(row_id, fields=True)

    if response["type"] != "row":
        raise ValueError(
            f"Expected `row_id` to resolve to a row, instead, it resolves to a `{response['type']}`"
        )

    # ask parent for column names
    parent_response = describe_row(response["parentId"])

    if parent_response["type"] != "database":
        raise ValueError(
            f"Expected parent of `{row_id}` to resolve to a database, instead, it resolves to a `{parent_response['type']}`"
        )

    # make a dictionary from column IDs to column names
    column_name_mapper = dict()
    column_cardinality_mapper = dict()
    for col in parent_response["cols"]:
        column_name_mapper[col["id"]] = col["name"]
        column_cardinality_mapper[col["id"]] = col["cardinality"]

    # now use this to construct the required dictionary
    row_data = dict()
    for field in response["fields"]:
        column_id = field["columnId"]
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
