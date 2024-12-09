import functools
import os
import pathlib
from typing import Optional

import confuse
from deeporigin.utils.core import _ensure_do_folder, in_aws_lambda

from ..exceptions import DeepOriginException

CONFIG_DIR = pathlib.Path(__file__).parent
DEFAULT_CONFIG_FILENAME = os.path.join(CONFIG_DIR, "default.yml")
CONFIG_YML_LOCATION = _ensure_do_folder() / "config.yml"

__all__ = ["get_value"]


# validate configuration
TEMPLATE = {
    "organization_id": confuse.String(),
    "bench_id": confuse.Optional(confuse.String()),
    "env": confuse.String(),
    "api_endpoint": confuse.Optional(confuse.String()),
    "nucleus_api_route": confuse.String(),
    "auth_domain": confuse.String(),
    "auth_device_code_endpoint": confuse.String(),
    "auth_token_endpoint": confuse.String(),
    "auth_audience": confuse.String(),
    "auth_grant_type": confuse.String(),
    "auth_client_id": confuse.String(),
    "auth_client_secret": confuse.String(),
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


@functools.cache
def get_value(
    config_file_location: Optional[str] = None, override_values: tuple = None
) -> confuse.templates.AttrDict:
    """Get the configuration for the Deep Origin CLI

    Args:
        user_config_filename: path to the user's configuration file
        override_values: values to use to override the default value

    Returns:
        :obj:`confuse.templates.AttrDict`: configuration for the Deep Origin CLI
    """

    # if we're running on AWS lambda, read config values
    # from env and return those.
    if in_aws_lambda():
        return confuse.AttrDict(
            dict(
                organization_id=os.environ.get("DEEP_ORIGIN_ORGANIZATION_ID"),
                api_endpoint=os.environ.get("DEEP_ORIGIN_API_ENDPOINT"),
                nucleus_api_route=os.environ.get("DEEP_ORIGIN_NUCLEUS_API_ROUTE"),
                env=os.environ.get("DEEP_ORIGIN_ENV"),
            )
        )

    # read the default configuration
    value = confuse.Configuration("deep_origin", __name__)

    # read the default configuration
    value.set_file(DEFAULT_CONFIG_FILENAME, base_for_paths=True)

    # read configuration overrides from the user
    if config_file_location is None:
        config_file_location = CONFIG_YML_LOCATION
    if os.path.isfile(config_file_location):
        value.set_file(
            config_file_location,
            base_for_paths=True,
        )

    # read configuration from environment variables
    value.set_env(sep="__")

    # set overriding values
    if override_values is not None:
        for k, v in override_values:
            value.set({k: v})

    try:
        validated_value = value.get(TEMPLATE)
    except confuse.exceptions.ConfigTypeError as exception:
        detail = str(exception).replace("\n", "\n  ")
        key = detail.split(":")[0].strip()
        raise DeepOriginException(
            title="Invalid configuration",
            message=f"The Deep Origin CLI and Python client requires a valid configuration. The {key} field is not valid:\n {detail}",
            fix=f"To fix this issue, run `deeporigin config set {key} <value>`",
        )

    return validated_value
