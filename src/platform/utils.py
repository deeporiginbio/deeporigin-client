import inspect
from urllib.parse import urljoin

import do_sdk_platform
from beartype import beartype
from box import Box
from deeporigin.auth import get_tokens
from deeporigin.config import get_value
from deeporigin.utils.core import PersistentDict, _get_method, _redact_responses


@beartype
def _get_api_client(*, api_name: str, configure: bool = True):
    """return a configured client for the API we want to access

    Args:
        api_name (str): name of the API

    Returns:
        configured client
    """

    if configure:
        configuration = do_sdk_platform.configuration.Configuration(
            host=urljoin(get_value()["api_endpoint"], "/api"),
            access_token=get_tokens()["access"],
        )

        client = do_sdk_platform.ApiClient(configuration=configuration)
    else:
        client = do_sdk_platform.ApiClient()

    api_class = getattr(do_sdk_platform, api_name)
    client = api_class(api_client=client)
    return client


@beartype
def _get_client_methods(client) -> set:
    """utility function to get methods from the client that return raw responses from the server"""
    methods = set(
        [
            attr
            for attr in dir(client)
            if callable(getattr(client, attr))
            and not attr.startswith("_")
            and "without_preload_content" in attr
        ]
    )

    return methods


def _create_function(*, method_path: str, api_name: str):
    """utility function the dynamically creates functions
    that wrap low-level functions in the DeepOrigin data API"""

    # we're constructing a client solely for the purposes
    # of inspecting its methods and extracting
    # function signatures. So we don't need any
    # authentication

    client = _get_api_client(
        configure=False,
        api_name=api_name,
    )

    method = _get_method(client, method_path)

    signature = inspect.signature(method)

    def dynamic_function(
        *,
        client=None,
        **kwargs,
    ):
        if client is None:
            client = _get_api_client(api_name=api_name)
        method = _get_method(client, method_path)

        # call the low level API method
        response = method(**kwargs)

        if not isinstance(response, dict):
            response = response.json()

        if "data" in response.keys():
            response = response["data"]
            if isinstance(response, list):
                response = [Box(item) for item in response]
            else:
                response = Box(response)
        else:
            response = Box(response)

        return response

    # attach the signature of the underlying method to the
    # function so that IDEs can display it properly
    dynamic_function.__signature__ = signature

    return dynamic_function