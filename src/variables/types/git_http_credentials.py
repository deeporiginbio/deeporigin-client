import os
import typing
import urllib.parse
import warnings

import pydantic
import validators

from ...utils import expand_user
from ...warnings import DeepOriginWarning
from ..base_type import Variable

T = typing.TypeVar("T")


class GitHttpCredentials(Variable):
    """Git HTTP credentials"""

    class Meta:
        platform_id = "Git HTTP credentials"

    repository: str = pydantic.Field(validate_default=True)
    username: str = pydantic.Field(min_length=1, validate_default=True)
    password: str = pydantic.Field(min_length=0, validate_default=True)

    @pydantic.field_validator("repository")
    @classmethod
    def value_must_be_a_valid_url(cls, value: str) -> str:
        """Verify whether a value is a valid URL

        Args:
            value (:obj:`str`): value

        Returns:
            :obj:`str`: value

        Raises:
            :obj:`ValueError`: if the value is not a valid URL
        """
        if not isinstance(value, str):
            raise ValueError(f"Repository must be a URL. {value} is not a URL.")

        if not validators.url(value) and not validators.url(f"https://{value}"):
            raise ValueError(f"Repository must be a URL. {value} is not a URL.")

        return value

    def parse_repository(self) -> tuple[str, str]:
        """Get the net location and scheme for the repository

        Returns:
            :obj:`tuple`: net location and scheme for the repository
        """
        if validators.url(self.repository):
            repository_obj = urllib.parse.urlparse(self.repository)
        else:
            repository_obj = urllib.parse.urlparse(f"https://{self.repository}")

        repository_scheme = repository_obj.scheme

        repository = f"{repository_obj.netloc}{repository_obj.path}"
        if repository_obj.params:
            repository += ";" + repository_obj.params
        if repository_obj.query:
            repository += "?" + repository_obj.query
        if repository_obj.fragment:
            repository += "#" + repository_obj.fragment

        return repository, repository_scheme

    def install(
        self,
        prev_self: Variable,
        user_home_dirname: str = os.path.expanduser("~"),
        overwrite: bool = False,
    ) -> None:
        """Install a variable/secret: Append the credentials to `~/.git-credentials`

        Args:
            prev_self (:obj:`Variable`): previous value of the variable
            user_home_dirname (:obj:`str`, optional): user's home directory
            overwrite (:obj:`bool`, optional): Whether to overwrite direct changes to variables and secrets
        """
        other_config, do_config, _ = self.read_git_credentials(self, user_home_dirname)

        prev_match = None
        if prev_self:
            prev_repository, prev_repository_scheme = prev_self.parse_repository()
            for credential in do_config:
                if credential.startswith(
                    f"{prev_repository_scheme}://{prev_self.username}"
                ) and credential.endswith(f"@{prev_repository}\n"):
                    prev_credential = f"{prev_repository_scheme}://{prev_self.username}:{prev_self.password}@{prev_repository}\n"
                    if credential == prev_credential:
                        prev_match = True
                    else:
                        prev_match = False

                    if overwrite or prev_match:
                        do_config.remove(credential)

        repository, repository_scheme = self.parse_repository()
        new_do_config = (
            f"{repository_scheme}://{self.username}:{self.password}@{repository}\n"
        )
        if overwrite or prev_self is None or prev_match:
            if new_do_config not in do_config:
                do_config.append(new_do_config)
        else:
            warnings.warn(
                (
                    f"Variable {self.name or ''} was not uninstalled because it has been modified directly. "
                    "To update its value, rerun with the `overwrite` option."
                ),
                DeepOriginWarning,
            )

        self.write_git_credentials(other_config, do_config, user_home_dirname)

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
        repository, repository_scheme = self.parse_repository()
        other_config, do_config, _ = self.read_git_credentials(self, user_home_dirname)

        uninstalled = False

        for credential in do_config:
            if (
                overwrite
                and credential.startswith(f"{repository_scheme}://{self.username}")
                and credential.endswith(f"@{repository}\n")
            ) or (
                credential
                == f"{repository_scheme}://{self.username}:{self.password}@{repository}\n"
            ):
                uninstalled = True
                do_config.remove(credential)

        if uninstalled:
            self.write_git_credentials(other_config, do_config, user_home_dirname)

        return uninstalled

    @classmethod
    def read_git_credentials(
        cls: T,
        variable: typing.Optional[T] = None,
        user_home_dirname: str = os.path.expanduser("~"),
    ) -> tuple[list[str], list[str], bool]:
        """Read a Git credentials file

        Args:
            variable (:obj:`GitHttpCredentials`, optional): variable
            user_home_dirname (:obj:`str`, optional): user's home directory

        Returns:
            :obj:`tuple`: regular Git credentials, Deep Origin Git credentials, whether the variable has already been installed
        """
        credentials_filename = expand_user(
            os.path.join("~", ".git-credentials"), user_home_dirname
        )
        os.makedirs(os.path.dirname(credentials_filename), exist_ok=True)

        if variable:
            repository, repository_scheme = variable.parse_repository()
            username = variable.username
            password = variable.password
            new_do_config = (
                f"{repository_scheme}://{username}:{password}@{repository}\n"
            )

        other_config = []
        do_config = []
        in_do_config = False
        if os.path.isfile(credentials_filename):
            with open(credentials_filename, "r") as file:
                in_deep_origin_cli_config = 0
                for line in file:
                    if line == "######## BEGIN DEEP ORIGIN CLI ########\n":
                        in_deep_origin_cli_config += 1
                    elif line == "######## END DEEP ORIGIN CLI ########\n":
                        in_deep_origin_cli_config -= 1
                    elif in_deep_origin_cli_config == 0:
                        other_config.append(line)
                    else:
                        if variable and line == new_do_config:
                            in_do_config = True
                        do_config.append(line)

        return (other_config, do_config, in_do_config)

    @classmethod
    def write_git_credentials(
        cls,
        other_config: list[str],
        do_config: list[str],
        user_home_dirname: str = os.path.expanduser("~"),
    ) -> None:
        """Write Git credentials to a file

        Args:
            other_config (:obj:`list`, optional): regular Git credentials
            do_config (:obj:`list`, optional): Deep Origin Git credentials
            user_home_dirname (:obj:`str`, optional): user's home directory
        """
        credentials_filename = expand_user(
            os.path.join("~", ".git-credentials"), user_home_dirname
        )
        os.makedirs(os.path.dirname(credentials_filename), exist_ok=True)

        with open(credentials_filename, "w") as file:
            file.write("".join(other_config).rstrip())

            file.write("\n\n")

            file.write("######## BEGIN DEEP ORIGIN CLI ########\n")
            file.write("".join(do_config).strip() + "\n")
            file.write("######## END DEEP ORIGIN CLI ########\n")
