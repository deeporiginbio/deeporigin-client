import json
import os
from dataclasses import dataclass
from typing import Literal, Union
from urllib.parse import parse_qs, urljoin, urlparse

import requests
from beartype import beartype
from deeporigin.config import get_value
from deeporigin.exceptions import DeepOriginException
from tabulate import tabulate

__all__ = [
    "expand_user",
]


RowType = Literal["row", "database", "workspace"]
"""Type of a row"""

FileStatus = Literal["ready", "archived"]
"""Status of a file"""

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
"""Type of a column"""

DATAFRAME_ATTRIBUTE_KEYS = {
    "file_ids",
    "id",
    "reference_ids",
}


Cardinality = Literal["one", "many"]

IDFormat = Literal["human-id", "system-id"]
"""Format of an ID"""

DatabaseReturnType = Literal["dataframe", "dict"]
"""Return type of a database"""


@dataclass
class PREFIXES:
    """Prefixes for CLI and Python client"""

    DO = "do://"
    FILE = "_file"
    DB = "_database"
    ROW = "_row"
    FOLDER = "_workspace"


@beartype
def _print_tree(tree: dict, offset: int = 0) -> None:
    """Helper function to pretty print a tree"""
    print(" " * offset + tree["hid"])

    if "children" not in tree.keys():
        return
    for child in tree["children"]:
        _print_tree(child, offset + 2)


@beartype
def _truncate(txt: str) -> str:
    """Utility function for truncating text"""

    TERMINAL_WIDTH, _ = os.get_terminal_size()
    txt = (
        (txt[: int(TERMINAL_WIDTH / 2)] + "…")
        if len(txt) > int(TERMINAL_WIDTH / 2)
        else txt
    )
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
            raise DeepOriginException(message=f"Failed to download file from {url}")

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
