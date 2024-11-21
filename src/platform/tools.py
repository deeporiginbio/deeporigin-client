"""bridge module to interact with the platform tools api"""

import sys

from deeporigin.platform.utils import add_functions_to_module

methods = add_functions_to_module(
    module=sys.modules[__name__],
    api_name="ToolsApi",
)

__all__ = list(methods)
