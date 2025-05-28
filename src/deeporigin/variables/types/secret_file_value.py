"""class to handle secret files"""

import abc
import os
import typing

import pydantic

from ..base_type import Variable
from .file import File

T = typing.TypeVar("T")


class SecretFileValue(Variable):
    """A value of a file"""

    @classmethod
    @property
    @abc.abstractmethod
    def FILENAME(self) -> str:
        pass  # pragma: no cover

    value: str = pydantic.Field(default=None, min_length=0, validate_default=True)

    def install(
        self,
        prev_self: Variable,
        user_home_dirname: str = os.path.expanduser("~"),
        overwrite: bool = False,
    ) -> None:
        """Install a variable/secret: Save contents to the specified location

        Args:
            prev_self (:obj:`Variable`): previous value of the variable
            user_home_dirname (:obj:`str`, optional): user's home directory
            overwrite (:obj:`bool`, optional): Whether to overwrite direct changes to variables and secrets
        """
        File.install_file(
            self.name,
            self.FILENAME,
            self.value,
            prev_self.value if prev_self else None,
            user_home_dirname=user_home_dirname,
            overwrite=overwrite,
        )

    def uninstall(
        self: T,
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
        return File.uninstall_file(
            self.name,
            self.FILENAME,
            self.value,
            user_home_dirname=user_home_dirname,
            overwrite=overwrite,
        )
