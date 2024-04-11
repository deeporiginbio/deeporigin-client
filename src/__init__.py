import os
import pathlib

from .do_api import (  # noqa: F401
    cache_do_api_tokens,
    get_do_api_tokens,
    read_cached_do_api_tokens,
    sign_into_do_platform,
)

__all__ = [
    "__version__",
    "get_do_api_tokens",
    "sign_into_do_platform",
    "cache_do_api_tokens",
    "read_cached_do_api_tokens",
]

SRC_DIR = pathlib.Path(__file__).parent

version_filename = os.path.join(SRC_DIR, "VERSION")

with open(version_filename, "r") as file:
    __version__ = file.read()
