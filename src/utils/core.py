"""this module contains some core utility functions that in turn do not depend on anything else in this library"""

import base64
import hashlib
import json
import os
import shutil
from datetime import datetime
from pathlib import Path
from typing import Union

from beartype import beartype
from box import Box
from tabulate import tabulate


class PrettyDict(Box):
    """A dict subclass with a custom pretty-print representation."""

    def __repr__(self):
        """pretty print a dict"""
        return json.dumps(
            dict(self),
            indent=2,
            ensure_ascii=False,
        )

    def _repr_html_(self):
        """pretty print a dict"""
        self.__repr__()


def set_key_to_value(
    obj: dict,
    target_key: str,
    new_value,
):
    """utility function to set a key in a dict to a value

    Args:
        obj (dict): the dictionary
        target_key (str): the key to be set
        new_value: the value to be set

    """
    if isinstance(obj, dict):
        for key, value in obj.items():
            if key == target_key:
                obj[key] = new_value
            else:
                set_key_to_value(value, target_key, new_value)
    elif isinstance(obj, list):
        for item in obj:
            set_key_to_value(item, target_key, new_value)


@beartype
def hash_file(file_path: str | Path) -> str:
    """
    Hashes the contents of a file using sha256

    Args:
        file_path (Union[str, Path]): Path to the file to be hashed.

    Returns:
        str: The hexadecimal hash digest of the file's contents.
    """
    file_path = Path(file_path)
    hasher = hashlib.new("sha256")

    # Read and update hash in chunks to handle large files efficiently
    with file_path.open("rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            hasher.update(chunk)

    return hasher.hexdigest()


@beartype
def hash_strings(strings: list[str]) -> str:
    """
    Computes a SHA-256 hash for a list of strings in an order-insensitive manner.

    The function sorts the input list, joins the sorted strings using a null-character delimiter,
    and returns the hexadecimal digest of the SHA-256 hash.

    Parameters:
        strings (list[str]): A list of strings.

    Returns:
        str: The hexadecimal SHA-256 hash of the sorted list.
    """
    sorted_strings = sorted(strings)

    # Use a delimiter
    combined_string = "--".join(sorted_strings)
    hash_obj = hashlib.sha256(combined_string.encode("utf-8"))
    return hash_obj.hexdigest()


@beartype
def _redact_responses(data: dict | list) -> dict | list:
    """utility function to redact responses of identifying information"""

    patterns = ["https://s3", "s3://data", "drn:identity::"]
    for pattern in patterns:
        data = _redact_data_in_dict(data, pattern)

    return data


def _redact_data_in_dict(
    data: dict | list,
    pattern: str = "https://s3",
):
    """Utility function to redacts patterns in a list or dict"""
    if isinstance(data, dict):
        return {
            key: _redact_data_in_dict(value, pattern) for key, value in data.items()
        }
    elif isinstance(data, list):
        return [_redact_data_in_dict(item, pattern) for item in data]
    elif isinstance(data, str) and pattern in data:
        return "<redacted>"

    else:
        # Leave other data types untouched

        return data


class PersistentDict:
    """
    A class that behaves like a dictionary but persists its data to a JSON file.
    """

    def __init__(self, file_path):
        """
        Initialize the PersistentDict.

        Args:
            file_path (str or Path): The path to the JSON file.
        """
        self.file_path = Path(file_path)
        self._data = self._load_or_initialize()

    def _load_or_initialize(self):
        """
        Load data from the JSON file if it exists, or initialize an empty dict.
        """
        if self.file_path.exists():
            with self.file_path.open("r", encoding="utf-8") as f:
                return json.load(f)
        else:
            self.file_path.parent.mkdir(
                parents=True, exist_ok=True
            )  # Ensure the directory exists
            self._save({})
            return {}

    def _save(self, data):
        """
        Save the current data to the JSON file.
        """
        with self.file_path.open("w", encoding="utf-8") as f:
            json.dump(data, f, indent=4)

    def __getitem__(self, key):
        return self._data[key]

    def __setitem__(self, key, value):
        self._data[key] = value
        self._save(self._data)

    def __delitem__(self, key):
        del self._data[key]
        self._save(self._data)

    def __contains__(self, key):
        return key in self._data

    def __iter__(self):
        return iter(self._data)

    def __len__(self):
        return len(self._data)

    def __repr__(self):
        return repr(self._data)

    def keys(self):
        """method to return the keys of this dict"""
        return self._data.keys()

    def values(self):
        """method to return the values of the underlying dict"""
        return self._data.values()

    def items(self):
        """method to return the items of the underlying dict"""
        return self._data.items()

    def get(self, key, default=None):
        """method to get a value by key"""
        return self._data.get(key, default)

    def update(self, *args, **kwargs):
        """method to update a value in the underlying dict"""
        self._data.update(*args, **kwargs)
        self._save(self._data)

    def clear(self):
        """method to clear the underlying dict"""
        self._data.clear()
        self._save(self._data)


def _get_method(obj, method_path: str):
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
def find_last_updated_row(rows: list[dict]) -> dict | None:
    """utility function to find the most recently updated row and return that object"""

    if len(rows) == 0:
        return None

    most_recent_date = None
    most_recent_row = rows[0]

    # Iterate over the list of objects
    for row in rows:
        current_date = datetime.strptime(row.dateUpdated, "%Y-%m-%d %H:%M:%S.%f")

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
