import collections.abc
import functools
import os
import pathlib

import confuse

from ..exceptions import DeepOriginException

CONFIG_DIR = pathlib.Path(__file__).parent
CONFIG_YML_LOCATION = os.path.expanduser(
    os.path.join(
        "~",
        ".deeporigin",
        "config.yml",
    )
)

__all__ = ["get_value"]


@functools.cache
def get_value(
    user_config_filenames: collections.abc.Iterable[str] = (
        CONFIG_YML_LOCATION,
        os.path.join(".deeporigin", "config.yml"),
    ),
) -> confuse.templates.AttrDict:
    """Get the configuration for the Deep Origin CLI

    Args:
        user_config_filenames (:obj:`list`, optional): paths to the user's configuration files [default: :obj:`['~/.deeporigin/config.yml', './.deeporigin/config.yml']`]

    Returns:
        :obj:`confuse.templates.AttrDict`: configuration for the Deep Origin CLI
    """
    value = confuse.Configuration("deep_origin", __name__)

    # read the default configuration
    default_config_filename = os.path.join(CONFIG_DIR, "default.yml")
    value.set_file(default_config_filename, base_for_paths=True)

    # read configuration overrides from the user
    for user_config_filename in user_config_filenames:
        if os.path.isfile(user_config_filename):
            value.set_file(user_config_filename, base_for_paths=True)
            break

    # read configuration from environment variables
    value.set_env(sep="__")

    # validate configuration
    template = {
        "organization_id": confuse.String(),
        "bench_id": confuse.String(),
        "env": confuse.String(),
        "api_endpoint": confuse.Optional(confuse.String()),
        "nucleus_api_route": confuse.String(),
        "graphql_api_route": confuse.String(),
        "auth_domain": confuse.String(),
        "auth_device_code_endpoint": confuse.String(),
        "auth_token_endpoint": confuse.String(),
        "auth_audience": confuse.String(),
        "auth_grant_type": confuse.String(),
        "auth_client_id": confuse.String(),
        "auth_client_secret": confuse.String(),
        "list_bench_variables_query_template": confuse.String(),
        "api_tokens_filename": confuse.Filename(),
        "variables_cache_filename": confuse.Filename(),
        "auto_install_variables_filename": confuse.Filename(),
        "feature_flags": confuse.Optional(
            confuse.MappingTemplate(
                {
                    "variables": confuse.TypeTemplate(bool),
                }
            )
        ),
    }
    try:
        validated_value = value.get(template)
    except confuse.exceptions.ConfigTypeError as exception:
        detail = str(exception).replace("\n", "\n  ")
        raise DeepOriginException(f"The configuration is not valid:\n  {detail}")

    return validated_value
