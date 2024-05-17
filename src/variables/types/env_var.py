import collections
import os
import typing
import warnings

import dotenv
import pydantic

from ...utils import expand_user
from ...warnings import DeepOriginWarning
from ..base_type import Variable

T = typing.TypeVar("T")


class EnvironmentVariable(Variable):
    """Environment variable"""

    class Meta:
        platform_id = "Non-secret environment variables"

    key: str = pydantic.Field(
        pattern=r"^[_a-zA-Z][_a-zA-Z0-9]*$", validate_default=True
    )
    value: str = pydantic.Field(default=None, min_length=0, validate_default=True)

    def install(
        self,
        prev_self: Variable,
        user_home_dirname: str = os.path.expanduser("~"),
        overwrite: bool = False,
    ) -> None:
        """Install a variable/secret: Append export command to `~/.bashrc`

        Args:
            prev_self (:obj:`Variable`): previous value of the variable
            user_home_dirname (:obj:`str`, optional): user's home directory
            overwrite (:obj:`bool`, optional): Whether to overwrite direct changes to variables and secrets
        """
        self.install_env_var(
            self.name,
            self.key,
            self.value,
            prev_self.value if prev_self else None,
            is_user_variable=True,
            user_home_dirname=user_home_dirname,
            overwrite=overwrite,
        )

    @classmethod
    def install_env_var(
        cls,
        name: str,
        key: str,
        value: str,
        prev_var_value: typing.Optional[str] = None,
        is_user_variable: bool = True,
        user_home_dirname: str = os.path.expanduser("~"),
        overwrite: bool = False,
    ) -> None:
        """Install a variable/secret: Append export command to `~/.bashrc`

        Args:
            name (:obj:`str`): name of the variable
            filename (:obj:`str`): filename of the variable
            value (:obj:`str`): value of the variable
            prev_var_value (:obj:`str`, optional): previous value of the variable
            is_user_variable (:obj:`bool`, optional): whether the variable is a user or system variable
            user_home_dirname (:obj:`str`, optional): user's home directory
            overwrite (:obj:`bool`, optional): Whether to overwrite direct changes to variables and secrets
        """
        env_filename = cls.get_env_filename(
            is_user_variable=is_user_variable, user_home_dirname=user_home_dirname
        )

        # Save variable to a dotenv file
        env = cls.read_deep_origin_environment_variables(env_filename)

        if overwrite or prev_var_value is None:
            prev_value = None
        else:
            prev_value = env.get(key, None)

        if overwrite or prev_value is None or prev_value == prev_var_value:
            env[key] = value
            cls.write_deep_origin_environment_variables(env, env_filename)

            # Add variable to bash profile
            cls.install_deep_origin_environment_variables(user_home_dirname)
        else:
            warnings.warn(
                (
                    f"Variable {name or ''} was not modified because it has been modified directly. "
                    "To update its value, rerun with the `overwrite` option."
                ),
                DeepOriginWarning,
            )

    def uninstall(
        self,
        user_home_dirname: str = os.path.expanduser("~"),
        overwrite: bool = False,
    ) -> bool:
        """Uninstall a variable/secret

        Args:
            user_home_dirname (:obj:`str`, optional): user's home directory
            overwrite (:obj:`bool`, optional): Whether to overwrite direct changes to variables and secrets

        Returns:
            :obj:`bool`: Whether the variable was uninstalled
        """
        return self.uninstall_env_var(
            self.name,
            self.key,
            self.value,
            is_user_variable=True,
            user_home_dirname=user_home_dirname,
            overwrite=overwrite,
        )

    @classmethod
    def uninstall_env_var(
        cls,
        name: str,
        key: str,
        value: str,
        is_user_variable: bool = True,
        user_home_dirname: str = os.path.expanduser("~"),
        overwrite: bool = False,
    ) -> bool:
        """Uninstall a variable/secret

        Args:
            name (:obj:`str`): name of the variable
            filename (:obj:`str`): filename of the variable
            value (:obj:`str`): value of the variable
            is_user_variable (:obj:`bool`, optional): whether the variable is a user or system variable
            user_home_dirname (:obj:`str`, optional): user's home directory
            overwrite (:obj:`bool`, optional): Whether to overwrite direct changes to variables and secrets

        Returns:
            :obj:`bool`: Whether the variable was uninstalled
        """
        env_filename = cls.get_env_filename(
            is_user_variable=is_user_variable, user_home_dirname=user_home_dirname
        )

        # Save variable to a dotenv file
        env = cls.read_deep_origin_environment_variables(env_filename)

        uninstalled = False

        if overwrite or env.get(key, None) == value:
            if key in env:
                uninstalled = True

                env.pop(key)

                cls.write_deep_origin_environment_variables(env, env_filename)

                # Add variable to bash profile
                cls.install_deep_origin_environment_variables(user_home_dirname)
        else:
            warnings.warn(
                (
                    f"Variable {name or ''} was not uninstalled because it has been modified directly. "
                    "To update its value, rerun with the `overwrite` option."
                ),
                DeepOriginWarning,
            )

        return uninstalled

    @classmethod
    def read_deep_origin_environment_variables(
        cls,
        env_filename: str,
    ) -> dict[str, str]:
        """Read the Deep Origin environment variables from a .env file

        Args:
            env_filename (:obj:`str`, optional): environment file name

        Returns:
            :obj:`dict`: environment variables
        """
        os.makedirs(os.path.dirname(env_filename), exist_ok=True)
        if os.path.isfile(env_filename) and os.stat(env_filename).st_size != 0:
            env = dotenv.dotenv_values(env_filename)
        else:
            env = collections.OrderedDict()
        return env

    @classmethod
    def write_deep_origin_environment_variables(
        cls, env: dict[str, str], env_filename: str
    ) -> None:
        """Write the Deep Origin environment variables to a .env file

        Args:
            env (:obj:`dict`): environment variables
            env_filename (:obj:`str`, optional): environment file name
        """
        with open(env_filename, "w") as file:
            for key, value in env.items():
                escaped_value = value.replace('"', '\\"')
                file.write(f'{key}="{escaped_value}"\n')

    @classmethod
    def install_deep_origin_environment_variables(
        cls,
        user_home_dirname: str = os.path.expanduser("~"),
    ) -> None:
        """Add sourcing Deep Origin environment variable into the user's bash profile

        Args:
            user_home_dirname (:obj:`str`, optional): user's home directory
        """
        user_env_filename = cls.get_env_filename(
            is_user_variable=True, user_home_dirname=user_home_dirname
        )
        system_env_filename = cls.get_env_filename(
            is_user_variable=False, user_home_dirname=user_home_dirname
        )

        if not os.path.isfile(user_env_filename):
            with open(user_env_filename, "w"):
                pass
        if not os.path.isfile(system_env_filename):
            with open(system_env_filename, "w"):
                pass

        bashrc_filename = expand_user(os.path.join("~", ".bashrc"), user_home_dirname)
        os.makedirs(os.path.dirname(bashrc_filename), exist_ok=True)

        env_file_installed = False
        if os.path.isfile(bashrc_filename):
            with open(bashrc_filename, "r") as file:
                for line in file:
                    if line == "######## BEGIN DEEP ORIGIN CLI ########\n":
                        env_file_installed = True
                        break

        if not env_file_installed:
            with open(bashrc_filename, "a") as file:
                file.write("\n\n")
                file.write("######## BEGIN DEEP ORIGIN CLI ########\n")
                file.write("# Load variables from the Deep Origin platform\n")
                file.write("set -o allexport\n")
                file.write(f"source {user_env_filename}\n")
                file.write(f"source {system_env_filename}\n")
                file.write("set +o allexport\n")
                file.write("######## END DEEP ORIGIN CLI ########\n")

    @classmethod
    def get_env_filename(
        cls,
        is_user_variable: bool = True,
        user_home_dirname: str = os.path.expanduser("~"),
    ) -> str:
        """Get the filename where platform managed environment variables should be stored

        Args:
            is_user_variable (:obj:`bool`, optional): whether the variable is a user or system variable
            user_home_dirname (:obj:`str`, optional): user's home directory

        Returns:
            :obj:`str`: Environment variable filename
        """
        if is_user_variable:
            return os.path.join(user_home_dirname, ".deeporigin", "user-variables.env")
        else:
            return os.path.join(
                user_home_dirname, ".deeporigin", "system-variables.env"
            )
