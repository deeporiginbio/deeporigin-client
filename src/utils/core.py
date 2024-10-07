"""this module contains some core utility functions that in turn do not depend on anything else in this library"""

import base64
import json
import os
import shutil
import typing
from datetime import datetime
from pathlib import Path

from beartype import beartype
from beartype.typing import List, Union
from tabulate import tabulate

T = typing.TypeVar("T")


def _get_method(obj, method_path):
    # Split the method path into components
    methods = method_path.split(".")

    # Traverse the attributes to get to the final method
    for method in methods:
        obj = getattr(obj, method)

    return obj


def humanize_file_size(file_size: int):
    """humanizes a file size in bytes"""
    for unit_prefix in ["", "K", "M", "G", "T", "P", "E", "Z", "Y"]:
        if file_size < 1024.0:
            return f"{file_size:.2f} {unit_prefix}B"
        file_size /= 1024.0


def sha256_checksum(file_path):
    """compute sha256 hash of a file"""

    import hashlib

    sha256_hash = hashlib.sha256()
    with open(file_path, "rb") as f:
        for byte_block in iter(lambda: f.read(4096), b""):
            sha256_hash.update(byte_block)
    return base64.b64encode(sha256_hash.digest()).decode()


def in_aws_lambda():
    """Returns True if running in AWS Lambda"""

    return "AWS_LAMBDA_FUNCTION_NAME" in os.environ


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
def _ensure_do_folder() -> Path:
    """makes sure that the deeporigin scratch folder exists

    If it doesn't exist, the folder is created. The location of this folder is typically in the home directory, but may be elsewhere if running in AWS Lambda"""

    if in_aws_lambda():
        deeporigin_path = Path("/tmp/.deeporigin")
    else:
        deeporigin_path = Path.home() / ".deeporigin"

    if not deeporigin_path.exists():
        deeporigin_path.mkdir(parents=True)

    return deeporigin_path


@beartype
def _get_api_tokens_filepath() -> Path:
    """get location of the api tokens file"""

    return _ensure_do_folder() / "api_tokens"


@beartype
def read_cached_tokens() -> dict:
    """Read cached API tokens"""

    with open(_get_api_tokens_filepath(), "r") as file:
        tokens = json.load(file)
    return tokens
