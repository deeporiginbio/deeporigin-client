"""bridge module to interact with the platform tools api"""

import sys

from deeporigin.platform.utils import _add_functions_to_module

__all__ = _add_functions_to_module(
    module=sys.modules[__name__],
    api_name="FilesApi",
    sdk_name="files",
)
