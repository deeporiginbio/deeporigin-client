"""bridge module to interact with the platform tools api"""

import os
import sys

from beartype import beartype

from deeporigin.platform.utils import _add_functions_to_module

__all__ = _add_functions_to_module(
    module=sys.modules[__name__],
    api_name="FilesApi",
    sdk_name="files",
)


@beartype
def upload_file(*, local_path: str, remote_path: str):
    """upload a file using the platform files api"""

    with open(local_path, "rb") as f:
        file_bytes = f.read()

    file_name = os.path.basename(local_path)

    # this function will be imported by _add_functions_to_module
    put_object(file_path=remote_path, file=(file_name, file_bytes))  # noqa: F821
