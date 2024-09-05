"""This module contains utility functions that are used internally by the Python Client and the CLI"""

import json
import os
import shutil
from dataclasses import dataclass
from datetime import datetime
from pathlib import Path
from typing import List, Literal, TypeVar, Union
from urllib.parse import parse_qs, urljoin, urlparse

import requests
from beartype import beartype
from deeporigin.config import get_value
from deeporigin.exceptions import DeepOriginException
from tabulate import tabulate

__all__ = [
    "expand_user",
]


T = TypeVar("T")

ObjectType = Literal["row", "database", "workspace"]
"""Type of a row. In Deep Origin, a row can be a database row, a database or a workspace"""

FileStatus = Literal["ready", "archived"]
"""Status of a file. Ready files are ready to be used, downloaded, and operated on."""

DataType = Literal[
    "integer",
    "str",
    "select",
    "date",
    "text",
    "file",
    "reference",
    "editor",
    "float",
    "boolean",
]
"""Type of a column in a Deep Origin database. See [this page in the documentation](https://docs.deeporigin.io/docs/os/data-hub/databases/columns) for more information."""

DATAFRAME_ATTRIBUTE_KEYS = {
    "file_ids",
    "id",
    "reference_ids",
}


Cardinality = Literal["one", "many"]
"""The cardinality defines whether a cell in a database can contain one or multiple objects"""

IDFormat = Literal["human-id", "system-id"]
"""Format of an ID"""

DatabaseReturnType = Literal["dataframe", "dict"]
"""Return type for [api.get_dataframe][src.data_hub.api.get_dataframe]"""


@dataclass
class PREFIXES:
    """Prefixes for CLI and Python client"""

    DO = "do://"
    FILE = "_file"
    DB = "_database"
    ROW = "_row"
    FOLDER = "_workspace"


def _get_pypi_version():
    """determines the latest version on PyPI"""

    response = requests.get("https://pypi.org/pypi/deeporigin/json")

    if response.status_code == 200:
        data = response.json()
        return data["info"]["version"]
    else:
        return None


@beartype
def construct_resource_url(
    *,
    name: str,
    row_type: ObjectType,
) -> str:
    """Constructs the URL for a resource

    Args:
        name (str): name of the resource
        row_type (ObjectType): type of the resource

    Returns:
        str: URL for the resource
    """

    env = get_value()["env"]
    org = get_value()["organization_id"]
    if env == "prod":
        url = f"https://os.deeporigin.io/org/{org}/data/{row_type}/{name}"
    else:
        url = f"https://os.{env}.deeporigin.io/org/{org}/data/{row_type}/{name}"

    return url


@beartype
def find_last_updated_row(rows: List[T]) -> T:
    """utility function to find the most recently updated row and return that object"""

    most_recent_date = None
    most_recent_row = rows[0]

    # Iterate over the list of objects
    for row in rows:
        current_date = datetime.strptime(row.date_updated, "%Y-%m-%d %H:%M:%S.%f")

        if most_recent_date is None or current_date > most_recent_date:
            most_recent_date = current_date
            most_recent_row = row

    return most_recent_row


@beartype
def _print_tree(tree: dict, offset: int = 0) -> None:
    """Helper function to pretty print a tree"""
    print(" " * offset + tree["hid"])

    if "children" not in tree.keys():
        return
    for child in tree["children"]:
        _print_tree(child, offset + 2)


def _truncate(txt: str) -> str:
    """Utility function for truncating text"""

    TERMINAL_WIDTH = int(shutil.get_terminal_size().columns / 2)

    if txt is None:
        return txt
    txt = str(txt)
    if len(txt) > TERMINAL_WIDTH:
        txt = txt[: TERMINAL_WIDTH - 3] + "..."
    return txt


@beartype
def _show_json(data: Union[list, dict]) -> None:
    """Utility for pretty printing JSON, used in the CLI"""

    print(json.dumps(data, indent=2))


@beartype
def _print_dict(
    data: dict,
    *,
    json: bool = True,
    transpose: bool = True,
    key_label: str = "Name",
) -> None:
    """Helper function to pretty print a dict as a table,
    used in the CLI"""

    if json:
        _show_json(data)
    else:
        if transpose:
            # truncate values so that long strings
            # don't break the table
            data = {key: _truncate(value) for key, value in data.items()}

            data = data.items()
            headers = [key_label, "Value"]
        else:
            headers = "keys"

        print(
            tabulate(
                data,
                headers=headers,
                tablefmt="rounded_outline",
            )
        )


@beartype
def _nucleus_url() -> str:
    """Returns URL for nucleus API endpoint"""
    url = urljoin(
        get_value()["api_endpoint"],
        get_value()["nucleus_api_route"],
    )
    if not url.endswith("/"):
        url += "/"

    return url


@beartype
def expand_user(path, user_home_dirname: str = os.path.expanduser("~")) -> str:
    """Expand paths that start with `~` by replacing it the user's home directory

    Args:
        path (:obj:`str`): path
        user_home_dirname (:obj:`str`, optional): user's home directory

    Returns:
        :obj:`str`: expanded path
    """
    if path == "~":
        return user_home_dirname
    elif path.startswith("~" + os.path.sep):
        return os.path.join(user_home_dirname, path[2:])
    else:
        return path


@beartype
def download_sync(url: str, save_path: str) -> None:
    """Concrete method to download a resource using GET and save to disk

    Args:
        url (str): url to download
        save_path (str): path to save file
    """

    with requests.get(url, stream=True) as response:
        if response.status_code != 200:
            raise DeepOriginException(
                message=f"File could not be downloaded from {url}"
            )

        with open(save_path, "wb") as file:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:  # Filter out keep-alive new chunks
                    file.write(chunk)


@beartype
def _parse_params_from_url(url: str) -> dict:
    """Utility function to extract params from a URL query

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


def _get_method(obj, method_path):
    # Split the method path into components
    methods = method_path.split(".")

    # Traverse the attributes to get to the final method
    for method in methods:
        obj = getattr(obj, method)

    return obj


def _ensure_do_folder():
    """makes sure that ~/.deeporigin exists"""

    deeporigin_path = Path.home() / ".deeporigin"

    if not deeporigin_path.exists():
        deeporigin_path.mkdir(parents=True)
