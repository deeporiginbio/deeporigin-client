import os
from urllib.parse import urljoin

from beartype import beartype
from deeporigin.config import get_value

__all__ = [
    "expand_user",
]


PREFIX = "deeporigin://"


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
