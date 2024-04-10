import abc
import os
import typing

import pydantic

from ..base_type import Variable
from .env_var import EnvironmentVariable

T = typing.TypeVar("T")


class SecretEnvironmentVariableValue(Variable):
    """A value of an environment variable"""

    @classmethod
    @property
    @abc.abstractmethod
    def KEY(self) -> str:
        pass  # pragma: no cover

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
        EnvironmentVariable.install_env_var(
            self.name,
            self.KEY,
            self.value,
            prev_self.value if prev_self else None,
            is_user_variable=False,
            user_home_dirname=user_home_dirname,
            overwrite=overwrite,
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
        return EnvironmentVariable.uninstall_env_var(
            self.name,
            self.KEY,
            self.value,
            is_user_variable=False,
            user_home_dirname=user_home_dirname,
            overwrite=overwrite,
        )
