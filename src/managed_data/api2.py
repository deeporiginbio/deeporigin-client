import inspect
import sys

from deeporigin import auth
from deeporigin._data_api import DeeporiginData

# this dict allows us to simply wrap endpoints in the low-level
# API so that we can more easily call them, without needing
# to explicitly write functions for them
WRAPPER_MAPPER = dict(
    create_file_download_url="create_file_download_url.create",
    describe_file="describe_file.retrieve",
    describe_row="describe_row.retrieve",
    list_database_rows="databases.rows.list",
    create_file_upload_url="create_file_upload.create",
    list_row_back_references="rows.back_references.list",
    list_mentions="mentions.list",
    describe_database_stats="describe_database_stats.retrieve",
)


def _get_default_client(client=None):
    """Internal function to instantiate client

    Creates and returns an authenticated client if
    not provided with one.

    Warning: Internal function
        Do not use this function

    Args:
        client: None, or a Client


    """
    if client is None:
        tokens = auth.get_tokens(refresh=False)
        access_token = tokens["access"]

        from deeporigin.config import get_value

        org_id = get_value()["organization_id"]

        client = DeeporiginData(
            bearer_token=access_token,
            org_id=org_id,
        )

    return client


def _get_method(obj, method_path):
    # Split the method path into components
    methods = method_path.split(".")

    # Traverse the attributes to get to the final method
    for method in methods:
        obj = getattr(obj, method)

    return obj


def _create_function(name, method_path):
    client = _get_default_client()
    method = _get_method(client, method_path)
    signature = inspect.signature(method)

    allowed_args = set(signature.parameters.keys())

    def dynamic_function(**kwargs):
        # Check for unexpected keyword arguments
        unexpected_args = set(kwargs.keys()) - allowed_args
        if unexpected_args:
            raise TypeError(
                f"{name}() got unexpected keyword arguments: {', '.join(unexpected_args)}"
            )

        # call the low level API method
        method(**kwargs)

    # attach the signature of the underlying method to the
    # function so that IDEs can display it properly
    dynamic_function.__signature__ = signature

    return dynamic_function


module_name = sys.modules[__name__]


# Dynamically create functions and attach them to the module
for name, method_path in WRAPPER_MAPPER.items():
    setattr(module_name, name, _create_function(name, method_path))
