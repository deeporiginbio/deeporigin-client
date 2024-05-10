import json
import os
from dataclasses import dataclass
from typing import Union
from urllib.parse import urljoin

from beartype import beartype
from deeporigin.config import get_value
from tabulate import tabulate

__all__ = [
    "expand_user",
]


@dataclass
class PREFIXES:
    """prefixes for CLI and python client"""

    DO = "do://"
    FILE = "_file"
    DB = "_database"
    ROW = "_row"
    WS = "_workspace"


@beartype
def _print_tree(tree: dict, offset: int = 0) -> None:
    """helper function to pretty print a tree"""
    print(" " * offset + tree["hid"])

    if "children" not in tree.keys():
        return
    for child in tree["children"]:
        _print_tree(child, offset + 2)


@beartype
def _truncate(txt: str) -> str:
    """utility function for truncating text"""

    TERMINAL_WIDTH, _ = os.get_terminal_size()
    txt = (
        (txt[: int(TERMINAL_WIDTH / 2)] + "…")
        if len(txt) > int(TERMINAL_WIDTH / 2)
        else txt
    )
    return txt


@beartype
def _show_json(data: Union[list, dict]) -> None:
    """utility for pretty printing JSON, used in the CLI"""

    print(json.dumps(data, indent=2))


@beartype
def _print_dict(
    data: dict,
    *,
    json: bool = True,
    transpose: bool = True,
) -> None:
    """helper function to pretty print a dict as a table,
    used in the CLI"""

    if json:
        _show_json(data)
    else:
        if transpose:
            data = data.items()
            headers = ["Name", "Value"]
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
    """returns URL for nucleus API endpoint"""
    return urljoin(
        get_value()["api_endpoint"],
        get_value()["nucleus_api_route"],
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
