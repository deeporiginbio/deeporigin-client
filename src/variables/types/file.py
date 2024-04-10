import os
import typing
import warnings

import pydantic

from ...utils import expand_user
from ...warnings import DeepOriginWarning
from ..base_type import Variable

T = typing.TypeVar("T")


class File(Variable):
    """File"""

    class Meta:
        platform_id = "Configuration file"

    filename: str = pydantic.Field(validate_default=True)
    value: str = pydantic.Field(default=None, min_length=0, validate_default=True)

    @pydantic.field_validator("filename")
    @classmethod
    def value_must_be_a_valid_path(cls, value: str) -> str:
        """Verify whether a value is a valid path

        Args:
            value (:obj:`str`): value

        Returns:
            :obj:`str`: value

        Raises:
            :obj:`ValueError`: if the value is not a valid URL
        """
        if not isinstance(value, str):
            raise ValueError(
                (
                    f"Filename must be an absolute path or begin with ~/ to "
                    "specify that a filename with your user home directory. "
                    f"{value} is not an absolute path and does not begin with ~/."
                )
            )

        expanded_value = expand_user(value, os.path.expanduser("~"))
        if not os.path.isabs(expanded_value):
            raise ValueError(
                (
                    f"Filename must be an absolute path or begin with ~/ to "
                    "specify that a filename with your user home directory. "
                    f"{value} is not an absolute path and does not begin with ~/."
                )
            )
        return value

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
        self.install_file(
            self.name,
            self.filename,
            self.value,
            prev_self.value if prev_self else None,
            user_home_dirname=user_home_dirname,
            overwrite=overwrite,
        )

    @classmethod
    def install_file(
        cls,
        var_name: str,
        var_filename: str,
        var_value: str,
        prev_var_value: typing.Optional[str] = None,
        user_home_dirname: str = os.path.expanduser("~"),
        overwrite: bool = False,
    ) -> None:
        """Install a variable/secret: Save contents to the specified location

        Args:
            var_name (:obj:`str`): name of the variable
            var_filename (:obj:`str`): filename of the variable
            var_value (:obj:`str`): value of the variable
            prev_var_value (:obj:`str`, optional): previous value of the variable
            user_home_dirname (:obj:`str`, optional): user's home directory
            overwrite (:obj:`bool`, optional): Whether to overwrite direct changes to variables and secrets
        """
        filename = expand_user(var_filename, user_home_dirname)
        value = var_value
        os.makedirs(os.path.dirname(filename), exist_ok=True)

        if overwrite or prev_var_value is None:
            prev_value = None
        elif os.path.isfile(filename):
            with open(filename, "r") as file:
                prev_value = file.read()
        else:
            prev_value = None

        if overwrite or prev_var_value is None or prev_value == prev_var_value:
            with open(filename, "w") as file:
                file.write(value)
        else:
            warnings.warn(
                (
                    f"Variable {var_name or ''} was not modified because it has been modified directly. "
                    "To update its value, rerun with the `overwrite` option."
                ),
                DeepOriginWarning,
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
        return self.uninstall_file(
            self.name,
            self.filename,
            self.value,
            user_home_dirname=user_home_dirname,
            overwrite=overwrite,
        )

    @classmethod
    def uninstall_file(
        cls,
        var_name: str,
        var_filename: str,
        var_value: str,
        user_home_dirname: str = os.path.expanduser("~"),
        overwrite: bool = False,
    ) -> bool:
        """Uninstall a variable/secret

        Args:
            var_name (:obj:`str`): name of the variable
            var_filename (:obj:`str`): filename of the variable
            var_value (:obj:`str`): value of the variable
            user_home_dirname (:obj:`str`, optional): user's home directory
            overwrite (:obj:`bool`, optional): Whether to overwrite direct changes to variables and secrets

        Returns:
            :obj:`bool`: Whether the variable was uninstalled
        """
        filename = expand_user(var_filename, user_home_dirname)
        value = var_value
        uninstalled = False
        if os.path.isfile(filename):
            if overwrite:
                prev_value = None
            else:
                with open(filename, "r") as file:
                    prev_value = file.read()

            if overwrite or prev_value == value:
                uninstalled = True
                os.remove(filename)
            else:
                warnings.warn(
                    (
                        f"Variable {var_name or ''} was not uninstalled because it has been modified directly. "
                        "To update its value, rerun with the `overwrite` option."
                    ),
                    DeepOriginWarning,
                )

        return uninstalled
