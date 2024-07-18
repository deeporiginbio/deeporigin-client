import inspect
import sys

from deeporigin import auth
from deeporigin.utils import _get_method
from deeporigin_data import DeeporiginData

# this dict allows us to simply wrap endpoints in the low-level
# API so that we can more easily call them, without needing
# to explicitly write functions for them
IGNORED_METHODS = {
    "close",
    "get",
    "post",
    "put",
    "delete",
    "get_api_list",
    "is_closed",
    "platform_headers",
    "patch",
    "with_options",
}


OVERRIDDEN_METHODS = {"list_rows", "convert_id_format"}


def _get_client_methods():
    # the only reason we're creating this client is to
    # extract methods from it. So no need to
    # authenticate
    client = DeeporiginData(token="", org_id="")
    methods = set(
        [
            attr
            for attr in dir(client)
            if callable(getattr(client, attr)) and not attr.startswith("_")
        ]
    ).difference(IGNORED_METHODS)

    return methods


def _get_default_client(client=None, refresh: bool = True):
    """Internal function to instantiate client

    Creates and returns an authenticated client if
    not provided with one.

    Warning: Internal function
        Do not use this function

    Args:
        client: None, or a Client


    """
    if client is None:
        tokens = auth.get_tokens(refresh=refresh)
        access_token = tokens["access"]

        import httpx
        from deeporigin.config import get_value

        value = get_value()

        org_id = value["organization_id"]
        base_url = httpx.URL.join(
            value["api_endpoint"],
            value["nucleus_api_route"],
        )

        client = DeeporiginData(
            token=access_token,
            org_id=org_id,
            base_url=base_url,
        )

    return client


def _create_function(method_path):
    """utility function the dynamically creates functions
    that wrap low-level functions in the DeepOrigin data API"""

    # we're constructing a client solely for the purposes
    # of inspecting its methods and extracting
    # function signatures. So we don't need any
    # authentication
    client = DeeporiginData(token="", org_id="")
    method = _get_method(client, method_path)

    signature = inspect.signature(method)

    def dynamic_function(*, client=None, **kwargs):
        if client is None:
            client = _get_default_client()
        method = _get_method(client, method_path)
        # call the low level API method
        response = method(**kwargs)
        if hasattr(response, "data"):
            return response.data
        else:
            # this should never happen
            return response

    # attach the signature of the underlying method to the
    # function so that IDEs can display it properly
    dynamic_function.__signature__ = signature

    return dynamic_function


module_name = sys.modules[__name__]


methods = _get_client_methods()
for method in methods:
    setattr(module_name, method, _create_function(method))


__all__ = list((methods | {"_get_default_client"}).difference(OVERRIDDEN_METHODS))
