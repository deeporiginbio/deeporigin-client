"""utility functions to fetch variables and secrets"""

import copy
import enum
import getpass
import os
import stat
import sys
import typing
import warnings

import crontab
import pydantic
import yaml

from deeporigin.platform.api import get_variables_and_secrets
from deeporigin.utils.core import _get_api_tokens_filepath

from ..config import get_value as get_config
from ..exceptions import DeepOriginException
from ..feature_flags import FeatureNotAvailableWarning
from ..feature_flags import get_value as get_feature_flags
from .base_type import Variable
from .types import VariableType

__all__ = [
    "VariableStatus",
    "get_variable_types_by_values",
    "install_variables",
    "get_variables_from_do_platform",
    "is_variable_modified",
    "enable_variable_auto_updating",
    "disable_variable_auto_updating",
    "uninstall_variables",
]

T = typing.TypeVar("T")


class VariableStatus(enum.Enum):
    """Type of change of a variable/secret"""

    unmodified = "unmodified"
    added = "added"
    modified = "modified"
    deleted = "deleted"


def get_variable_types_by_values(values: list[str]) -> list[VariableType]:
    """Get a list of types of variables and secrets by their values

    Args:
        values (:obj:`list`): values of the variables and secrets

    Returns:
        :obj:`list`: variables and secrets

    Raises:
        :obj:`ValueError`: if one or more of the values is invalid
    """
    types = []
    invalid_values = set()
    for value in values:
        try:
            type = VariableType[value]
            if type not in types:
                types.append(type)
        except KeyError:
            invalid_values.add(value)

    if invalid_values:
        msg = "\n\n".join(
            [
                "The following types of variables and secrets are not valid:\n  {}".format(
                    "\n  ".join(sorted(invalid_values))
                ),
                "The following types of variables and secrets are supported:\n  {}".format(
                    "\n  ".join(sorted(VariableType.__members__.keys()))
                ),
            ]
        )
        raise DeepOriginException(message=msg)

    return types


def install_variables(
    user: bool = True,
    org: bool = True,
    types: typing.Iterable[VariableType] = tuple(VariableType.__members__.values()),
    user_home_dirname: str = os.path.expanduser("~"),
    overwrite: bool = False,
) -> dict[str, dict[str, typing.Union[str, VariableType]]]:
    """
    Retrieve the variables and secrets for your workstation from the Deep Origin OS and
    install them into the workstation. Includes your variables and secrets, as well as those of
    the parent organization of the workstation.

    Args:
        user (:obj:`bool`, optional): Whether to retrieve and install user variables and secrets
        org (:obj:`bool`, optional): Whether to retrieve and install organization variables and secrets
        types (:obj:`typing.Iterable`, optional): Types of variables and secrets to retrieve and
        user_home_dirname (:obj:`str`, optional): user's home directory
        overwrite (:obj:`bool`, optional): Whether to overwrite direct changes to variables and secrets

    Returns:
        :obj:`dict`: dictionary that maps the names of variables to whether they were modified or not
    """
    prev_variables = get_variables_from_local()

    # get variables from the DO OS
    variables = get_variables_from_do_platform(user=user, org=org, types=types)

    # determine which variables were added, modified, removed, or unmodified
    variable_modification = {}
    for variable in variables:
        modified = is_variable_modified(variable)
        variable_key = get_variable_local_key(variable)
        if modified:
            prev_variable_serialized = prev_variables.get(variable_key, None)
            if prev_variable_serialized:
                prev_variable_serialized.pop("type")
                prev_variable = variable.__class__(**prev_variable_serialized)
                uninstalled = prev_variable.uninstall(
                    user_home_dirname=user_home_dirname,
                    overwrite=overwrite,
                )

                variable_modification[variable_key] = {
                    "status": VariableStatus.modified,
                    "type": variable.__class__.__name__,
                    "name": variable.name,
                }
            else:
                prev_variable = None
                uninstalled = False
                variable_modification[variable_key] = {
                    "status": VariableStatus.added,
                    "type": variable.__class__.__name__,
                    "name": variable.name,
                }
        else:
            prev_variable = variable
            uninstalled = False
            variable_modification[variable_key] = {
                "status": VariableStatus.unmodified,
                "type": variable.__class__.__name__,
                "name": variable.name,
            }

        variable.install(
            prev_variable,
            user_home_dirname=user_home_dirname,
            overwrite=overwrite or uninstalled,
        )

    # determine which variables were deleted
    cur_variables = get_variables_from_local()

    for variable in variables:
        variable_key = get_variable_local_key(variable)
        if variable_key in prev_variables:
            prev_variables.pop(variable_key)

    for serialized_variable in prev_variables.values():
        type = serialized_variable.pop("type")
        cls = VariableType[type].value
        variable = cls(**serialized_variable)
        variable_key = get_variable_local_key(variable)
        if variable_key in cur_variables:
            cur_variables.pop(variable_key)

        variable.uninstall(user_home_dirname=user_home_dirname, overwrite=overwrite)

        variable_modification[variable_key] = {
            "status": VariableStatus.deleted,
            "type": variable.__class__.__name__,
            "name": variable.name,
        }

    export_variables_to_local(cur_variables)

    # return information about which variables were modified
    return variable_modification


def get_variables_from_do_platform(
    user: bool = True,
    org: bool = True,
    types: typing.Iterable[VariableType] = tuple(VariableType.__members__.values()),
) -> list[Variable]:
    """
    Retrieve the variables and secrets for your workstation from the Deep Origin OS.
    Includes your variables and secrets, as well as those of the parent organization
    of the workstation.

    Args:
        user (:obj:`bool`, optional): Whether to retrieve and install user variables and secrets
        org (:obj:`bool`, optional): Whether to retrieve and install organization variables and secrets
        types (:obj:`typing.Iterable`, optional): Types of variables and secrets to retrieve and install

    Returns:
        :obj:`list`: variables
    """
    if not get_feature_flags().variables:
        msg = "Updating variables is not yet available. For beta access, please contact support at support@deeporigin.com."
        warnings.warn(msg, FeatureNotAvailableWarning, stacklevel=2)
        return []

    classes = tuple([type.value for type in types])

    secrets = get_variables_and_secrets()

    # reshape secrets to match what we had in the old
    # graph QL response
    reshaped_secrets = []
    keys = ["key", "filename"]
    for secret in secrets:
        data = dict()
        data["name"] = secret["attributes"]["name"]
        data["type"] = secret["attributes"]["type"]
        for key in keys:
            if key in secret["attributes"].keys():
                data[key] = secret["attributes"][key]
        data["value"] = secret["attributes"]["value"]
        data["drn"] = "drn:variable:" + secret["id"]
        reshaped_secrets.append(data)

    # Parse GraphQL response
    variables = []
    invalid_variables = []
    for serialized_variable in reshaped_secrets:
        try:
            variable = deserialize_variable(serialized_variable)

            if isinstance(variable, classes):
                variables.append(variable)
        except (DeepOriginException, pydantic.ValidationError) as exception:
            invalid_variables.append(
                (
                    f"Variable {serialized_variable['name']} ({serialized_variable['drn']}) is invalid:\n"
                    "  " + str(exception).strip().replace("\n", "\n  ")
                )
            )

    if invalid_variables:
        msg = f"{len(invalid_variables)} variables are not valid.\n\n" + "\n\n".join(
            invalid_variables
        )
        raise DeepOriginException(message=msg)

    # return variables
    return variables


def deserialize_variable(serialized: dict) -> Variable:
    """Deserialize a variable from an GraphQL representation of a variable

    Args:
        serialized (:obj:`dict`): GraphQL representation of a variable

    Returns:
        :obj:`Variable`: variable
    """
    serialized = serialized.copy()

    type_platform_id = serialized.pop("type")

    for cls in VariableType.__members__.values():
        if cls.value.Meta.platform_id == type_platform_id:
            return cls.value.from_platform(serialized)

    raise DeepOriginException(
        message=f"{type} is not a valid type of variable or secret."
    )


def is_variable_modified(variable: Variable) -> bool:
    """Determine whether a variable was modified and if so store the current value

    Args:
        variable (:obj:`Variable`): variable

    Returns:
        :obj:`bool`: whether the variable was modified
    """
    # Read the previous variables from the local store
    variables = get_variables_from_local()

    # Determine if the variable was modified
    key = get_variable_local_key(variable)
    variable_dict = variable.__dict__
    variable_dict["type"] = variable.__class__.__name__

    old_variable_dict = variables.get(key, {})
    modified = old_variable_dict != variable_dict

    # If the variable was modified, update the local store
    if modified:
        variables[key] = variable_dict
        export_variables_to_local(variables)

    # Return whether the variable was modified
    return modified


def get_variables_from_local() -> dict:
    """Read variables from a local file

    Returns:
        :obj:`dict`: variables
    """
    config = get_config()

    if os.path.isfile(config.variables_cache_filename):
        with open(config.variables_cache_filename, "r") as file:
            return yaml.load(file, Loader=yaml.CLoader)
    else:
        return {}


def export_variables_to_local(variables: dict) -> None:
    """Export variables to a local file

    Args:
        variables (:obj:`dict`): variables
    """
    config = get_config()

    os.makedirs(os.path.dirname(config.variables_cache_filename), exist_ok=True)
    with open(config.variables_cache_filename, "w") as file:
        yaml.dump(variables, file, Dumper=yaml.CDumper)


def get_variable_local_key(variable: Variable) -> str:
    """Get the local key for the variable

    Args:
        variable (:obj:`Variable`): variable

    Returns:
        :obj:`str`: unique identifier for the variable
    """
    return variable.drn


def enable_variable_auto_updating(
    user: bool = True,
    org: bool = True,
    types: typing.Iterable[VariableType] = tuple(VariableType.__members__.values()),
    time_period_min: int = 30,
    overwrite: bool = False,
    cli: typing.Optional[str] = None,
) -> None:
    """
    Configure the variables and secrets for your workstation from the Deep Origin
    OS to be automatically installed once added or modified. Includes your
    variables and secrets, as well as those of the parent organization of the workstation.

    Args:
        user (:obj:`bool`, optional): Whether to retrieve and install user variables and secrets
        org (:obj:`bool`, optional): Whether to retrieve and install organization variables and secrets
        types (:obj:`typing.Iterable`, optional): Types of variables and secrets to retrieve and install
        time_period_min (:obj:`int`, optional): Time period for updating variables from the Deep Origin OS [default: 30, units: min]
        overwrite (:obj:`bool`, optional): Whether to overwrite direct changes to variables and secrets
    """
    # login
    if not get_feature_flags().variables:
        msg = "Updating variables is not yet available. For beta access, please contact support at support@deeporigin.com."
        warnings.warn(msg, FeatureNotAvailableWarning, stacklevel=2)
        return []

    # load the CronTab configuration
    with crontab.CronTab(user=True) as cron_tab:
        # get the id of the cron job
        cron_job_id = get_auto_install_variables_cronjob_id()

        # remove previous job
        for job in cron_tab.find_comment(comment=cron_job_id):
            if os.path.isfile(job.command):
                os.remove(job.command)
        cron_tab.remove_all(comment=cron_job_id)

        # assemble the command for the job
        config = get_config()
        escaped_organization_id = config.organization_id.replace("'", "\\'")
        escaped_workstation_id = config.bench_id.replace("'", "\\'")
        escaped_env = config.env.replace("'", "\\'")
        escaped_api_endpoint = config.api_endpoint.replace("'", "\\'")
        escaped_auth_domain = config.auth_domain.replace("'", "\\'")
        escaped_auth_device_code_endpoint = config.auth_device_code_endpoint.replace(
            "'", "\\'"
        )
        escaped_auth_token_endpoint = config.auth_token_endpoint.replace("'", "\\'")
        escaped_auth_audience = config.auth_audience.replace("'", "\\'")
        escaped_auth_grant_type = config.auth_grant_type.replace("'", "\\'")
        escaped_auth_client_id = config.auth_client_id.replace("'", "\\'")
        escaped_auth_client_secret = config.auth_client_secret.replace("'", "\\'")

        escaped_api_tokens_filename = str(_get_api_tokens_filepath()).replace(
            "'", "\\'"
        )
        escaped_variables_cache_filename = config.variables_cache_filename.replace(
            "'", "\\'"
        )
        escaped_auto_install_variables_filename = (
            config.auto_install_variables_filename.replace("'", "\\'")
        )

        commands = []
        username = getpass.getuser()
        user_home = os.path.expanduser("~").replace("'", "\\'")
        path = os.getenv("PATH").replace("'", "\\'")
        commands.append(f"export USER='{username}'")
        commands.append(f"export HOME='{user_home}'")
        commands.append(f"export PATH='{path}'")
        commands.append(
            f"export DEEP_ORIGIN_ORGANIZATION_ID='{escaped_organization_id}'"
        )
        commands.append(f"export DEEP_ORIGIN_BENCH_ID='{escaped_workstation_id}'")
        commands.append(f"export DEEP_ORIGIN_ENV='{escaped_env}'")
        commands.append(f"export DEEP_ORIGIN_API_ENDPOINT='{escaped_api_endpoint}'")
        commands.append(f"export DEEP_ORIGIN_AUTH_DOMAIN='{escaped_auth_domain}'")
        commands.append(
            f"export DEEP_ORIGIN_AUTH_DEVICE_CODE_ENDPOINT='{escaped_auth_device_code_endpoint}'"
        )
        commands.append(
            f"export DEEP_ORIGIN_AUTH_TOKEN_ENDPOINT='{escaped_auth_token_endpoint}'"
        )
        commands.append(f"export DEEP_ORIGIN_AUTH_AUDIENCE='{escaped_auth_audience}'")
        commands.append(
            f"export DEEP_ORIGIN_AUTH_GRANT_TYPE='{escaped_auth_grant_type}'"
        )
        commands.append(f"export DEEP_ORIGIN_AUTH_CLIENT_ID='{escaped_auth_client_id}'")
        commands.append(
            f"export DEEP_ORIGIN_AUTH_CLIENT_SECRET='{escaped_auth_client_secret}'"
        )

        commands.append(
            f"export DEEP_ORIGIN_API_TOKENS_FILENAME='{escaped_api_tokens_filename}'"
        )
        commands.append(
            f"export DEEP_ORIGIN_VARIABLES_CACHE_FILENAME='{escaped_variables_cache_filename}'"
        )
        commands.append(
            f"export DEEP_ORIGIN_AUTO_INSTALL_VARIABLES_FILENAME='{escaped_auto_install_variables_filename}'"
        )

        if not cli:
            cli = sys.argv[0]
            if not cli.endswith("/deeporigin") and not cli.endswith(
                "pytest/__main__.py"
            ):
                raise DeepOriginException(
                    message="Auto installation must be run from the Deep Origin CLI"
                )

        cli_command = [cli, "variables", "install"]

        if not user:
            cli_command.append("--no-user")

        if not org:
            cli_command.append("--no-org")

        cli_command.append("--type")
        cli_command.extend([type.name for type in types])

        if overwrite:
            cli_command.append("--overwrite")

        log_filename = os.path.expanduser(
            os.path.join("~", ".log", "deeporigin", "variables-auto-install.log")
        )
        os.makedirs(os.path.dirname(log_filename), exist_ok=True)

        cli_command = " ".join(cli_command)
        commands.append(f"{cli_command} >> {log_filename} 2>&1")

        job_filename = config.auto_install_variables_filename
        with open(job_filename, "w") as file:
            file.write("#!/usr/bin/env bash\n")
            file.write("\n".join(commands))
            file.write("\n")
        os.chmod(
            job_filename,
            stat.S_IRUSR | stat.S_IWUSR | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH,
        )

        # install cron job
        job = cron_tab.new(command=job_filename, comment=cron_job_id)
        job.minute.every(time_period_min)


def disable_variable_auto_updating() -> None:
    """Disable the auto updating of Deep Origin variables"""
    with crontab.CronTab(user=True) as cron_tab:
        cron_job_id = get_auto_install_variables_cronjob_id()
        for job in cron_tab.find_comment(comment=cron_job_id):
            job.enable(False)


def get_auto_install_variables_cronjob_id() -> str:
    """Get the id for the cron job for auto updating Deep Origin variables

    Returns:
        :obj:`str`: id of the cron job
    """
    return f"deep-origin-install-variables-{getpass.getuser()}"


def uninstall_variables(
    user: bool = True,
    org: bool = True,
    types: typing.Iterable[VariableType] = tuple(VariableType.__members__.values()),
    user_home_dirname: str = os.path.expanduser("~"),
    overwrite: bool = False,
) -> None:
    """
    Remove variables and secrets retrieved from the Deep Origin OS, including
    uninstalling variables and secrets.

    Args:
        user (:obj:`bool`, optional): Whether to retrieve and install user variables and secrets
        org (:obj:`bool`, optional): Whether to retrieve and install organization variables and secrets
        types (:obj:`typing.Iterable`, optional): Types of variables and secrets to retrieve and install
        user_home_dirname (:obj:`str`, optional): user's home directory
        overwrite (:obj:`bool`, optional): Whether to overwrite direct changes to variables and secrets
    """
    # Read the previous variables from the local store
    serialized_variables = get_variables_from_local()

    for serialized_variable in list(serialized_variables.values()):
        serialized_variable = copy.copy(serialized_variable)
        type = VariableType[serialized_variable.pop("type")]
        if type in types:
            cls = type.value
            variable = cls(**serialized_variable)
            variable.uninstall(user_home_dirname=user_home_dirname, overwrite=overwrite)
            serialized_variables.pop(get_variable_local_key(variable))

    export_variables_to_local(serialized_variables)
