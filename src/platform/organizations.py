"""bridge module to interact with the platform organizations api"""

import sys

from deeporigin.platform.utils import add_functions_to_module

API_NAME = "OrganizationsApi"


methods = add_functions_to_module(
    module=sys.modules[__name__],
    api_name=API_NAME,
)

__all__ = list(methods)
