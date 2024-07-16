import inspect
import sys
from typing import Optional, Union

from deeporigin import auth
from deeporigin._data_api import DeeporiginData
from deeporigin.utils import _get_method

# this dict allows us to simply wrap endpoints in the low-level
# API so that we can more easily call them, without needing
# to explicitly write functions for them
WRAPPER_MAPPER = dict(
    create_file_download_url=dict(
        method="create_file_download_url.create",
    ),
    describe_file=dict(
        method="describe_file.retrieve",
    ),
    describe_row=dict(
        method="describe_row.retrieve",
    ),
    list_database_rows=dict(method="databases.rows.list"),
    create_file_upload_url=dict(
        method="create_file_upload.create",
    ),
    list_row_back_references=dict(
        method="rows.back_references.list",
    ),
    list_mentions=dict(
        method="mentions.list",
    ),
    describe_database_stats=dict(
        method="describe_database_stats.retrieve",
    ),
    convert_id_format=dict(
        method="convert_id_format.convert",
    ),
    delete_databse_column=dict(
        method="delete_database_column.delete",
    ),
    delete_rows=dict(method="delete_rows.delete"),
    create_database=dict(method="databases.create"),
    create_workspace=dict(method="workspaces.create"),
    update_workspace=dict(method="update_workspace.run"),
    update_database=dict(method="update_database.run"),
    add_database_column=dict(method="add_database_column.add"),
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
        tokens = auth.get_tokens(refresh=True)
        access_token = tokens["access"]

        from deeporigin.config import get_value

        org_id = get_value()["organization_id"]

        client = DeeporiginData(
            bearer_token=access_token,
            org_id=org_id,
        )

    return client


def _create_function(name, data, client=None):
    """utility function the dynamically creates functions
    that wrap low-level functions in the DeepOrigin data API"""

    method_path = data["method"]
    if client is None:
        client = _get_default_client()
    method = _get_method(client, method_path)

    signature = inspect.signature(method)

    def dynamic_function(*, client=None, **kwargs):
        if client is None:
            client = _get_default_client()
        method = _get_method(client, method_path)
        # call the low level API method
        return method(**kwargs)

    # attach the signature of the underlying method to the
    # function so that IDEs can display it properly
    dynamic_function.__signature__ = signature

    return dynamic_function


module_name = sys.modules[__name__]


def _generate_functions(client=None):
    # Dynamically create functions and attach them to the module
    for name, data in WRAPPER_MAPPER.items():
        setattr(module_name, name, _create_function(name, data))


# generate functions for the default client
_generate_functions()

__all__ = [function for function in WRAPPER_MAPPER.keys()]
