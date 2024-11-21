"""bridge module to interact with the platform organizations api"""

import sys

from deeporigin.platform.utils import (
    _create_function,
    _get_api_client,
    _get_client_methods,
)

API_NAME = "OrganizationsApi"


this_module_name = sys.modules[__name__]


methods = _get_client_methods(_get_api_client(api_name=API_NAME))


for method in methods:
    # clean up the name so that it's more readable
    sanitized_method_name = method.split("controller")[0].rstrip("_")

    # add this function as an attribute to this module
    # so that we can call it
    setattr(
        this_module_name,
        sanitized_method_name,
        _create_function(
            method_path=method,
            api_name=API_NAME,
        ),
    )


__all__ = list(methods)
