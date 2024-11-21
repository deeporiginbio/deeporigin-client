"""bridge module to interact with the platform organizations api"""

import sys

from deeporigin.platform.utils import add_functions_to_module

methods = add_functions_to_module(
    module=sys.modules[__name__],
    api_name="OrganizationsApi",
)

__all__ = list(methods)
