import abc
import os
import typing

import pydantic

T = typing.TypeVar("T")


class Variable(abc.ABC, pydantic.BaseModel):
    """Variable or secret"""

    class Meta:
        platform_id = "SSH private keys"

    drn: str = pydantic.Field(pattern=r"^drn:", validate_default=True)
    name: str = pydantic.Field(default="", min_length=1, validate_default=True)

    @classmethod
    def from_platform(cls: type[T], platform_value: dict) -> T:
        """Create an instance from a dictionary from the Deep Origin platform

        Args:
            platform_value (:obj:`dict`): platform representation of the variable/secret

        Returns:
            :obj:`type[T]`: Python instance of the variable/secret
        """
        return cls(**platform_value)

    @abc.abstractmethod
    def install(
        self: T,
        prev_self: T,
        user_home_dirname: str = os.path.expanduser("~"),
        overwrite: bool = False,
    ) -> None:
        """Install a variable/secret

        Args:
            prev_self (:obj:`Variable`): previous value of the variable
            user_home_dirname (:obj:`str`, optional): user's home directory
            overwrite (:obj:`bool`, optional): Whether to overwrite direct changes to variables and secrets
        """
        pass  # pragma: no cover

    @abc.abstractmethod
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
        pass  # pragma: no cover
