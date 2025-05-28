"""class to handle SSH keys"""

import os
import subprocess
import tempfile
import typing

import pydantic

from deeporigin.utils.core import expand_user

from .secret_file import SecretFile

T = typing.TypeVar("T")


class PrivateSshKey(SecretFile):
    """Private SSH key"""

    class Meta:
        platform_id = "SSH private keys"

    @pydantic.field_validator("value")
    @classmethod
    def value_must_be_a_ssh_private_key(cls, value: str) -> str:
        """Verify whether a value is a valid SSH private key

        Args:
            value (:obj:`str`): value

        Returns:
            :obj:`str`: value

        Raises:
            :obj:`ValueError`: if the value is not a valid SSH private key
        """
        if not isinstance(value, str):
            raise ValueError(
                (
                    f"Value must be an SSH private key. "
                    f"{value} is not an SSH private key."
                )
            )

        if not (
            value.startswith("-----BEGIN OPENSSH PRIVATE KEY-----")
            or value.startswith("-----BEGIN RSA PRIVATE KEY-----")
            or value.startswith("-----BEGIN DSA PRIVATE KEY-----")
            or value.startswith("-----BEGIN EC PRIVATE KEY-----")
            or value.startswith("-----BEGIN PRIVATE KEY-----")
            or value.startswith("-----BEGIN ENCRYPTED PRIVATE KEY-----")
        ):
            raise ValueError(
                (
                    f"Value must be an SSH private key. "
                    f"{value} is not an SSH private key."
                )
            )

        _, key_filename = tempfile.mkstemp()
        with open(key_filename, "w") as file:
            file.write(value)

        cmd = ["ssh-keygen", "-y", "-e", "-f", key_filename]

        completed_process = subprocess.run(
            cmd,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )

        os.remove(key_filename)

        if completed_process.returncode != 0:
            raise ValueError(
                (
                    f"Value must be an SSH private key. "
                    f"{value} is not an SSH private key."
                )
            )

        return value

    def install(
        self,
        prev_self: T,
        user_home_dirname: str = os.path.expanduser("~"),
        overwrite: bool = False,
    ) -> None:
        """Install a variable/secret: Save the key to the specified location, and declare in `~/.ssh/config`

        Args:
            prev_self (:obj:`Variable`): previous value of the variable
            user_home_dirname (:obj:`str`, optional): user's home directory
            overwrite (:obj:`bool`, optional): Whether to overwrite direct changes to variables and secrets
        """
        super().install(prev_self, user_home_dirname, overwrite)

        if overwrite or prev_self is None or self.value == prev_self.value:
            ssh_config_filename = expand_user(
                os.path.join("~", ".ssh", "config"), user_home_dirname
            )
            os.makedirs(os.path.dirname(ssh_config_filename), exist_ok=True)

            other_config, do_config, in_do_config = self.read_ssh_config(
                self, user_home_dirname
            )

            if not in_do_config:
                new_do_config = f"IdentityFile {self.filename}\n"
                do_config.append(new_do_config)

            self.write_ssh_config(other_config, do_config, user_home_dirname)

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
        uninstalled = super().uninstall(
            user_home_dirname=user_home_dirname, overwrite=overwrite
        )

        if uninstalled:
            other_config, do_config, in_do_config = self.read_ssh_config(
                self, user_home_dirname
            )

            if in_do_config:
                new_do_config = f"IdentityFile {self.filename}\n"
                do_config.remove(new_do_config)

            self.write_ssh_config(other_config, do_config, user_home_dirname)

        return uninstalled

    @classmethod
    def read_ssh_config(
        cls: T,
        variable: typing.Optional[T] = None,
        user_home_dirname: str = os.path.expanduser("~"),
    ) -> tuple[list[str], list[str], bool]:
        """Read an SSH configuration file

        Args:
            variable (:obj:`PrivateSshKey`, optional): variable
            user_home_dirname (:obj:`str`, optional): user's home directory

        Returns:
            :obj:`tuple`: regular SSH configuration, Deep Origin SSH configuration, whether the variable has already been installed
        """
        ssh_config_filename = expand_user(
            os.path.join("~", ".ssh", "config"), user_home_dirname
        )

        other_config = []
        do_config = []
        if variable:
            new_do_config = f"IdentityFile {variable.filename}\n"
        in_do_config = False
        if os.path.isfile(ssh_config_filename):
            with open(ssh_config_filename, "r") as file:
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
    def write_ssh_config(
        cls,
        other_config: list[str],
        do_config: list[str],
        user_home_dirname: str = os.path.expanduser("~"),
    ) -> None:
        """Write SSH configuration to a file

        Args:
            other_config (:obj:`list`, optional): regular SSH configuration
            do_config (:obj:`list`, optional): Deep Origin SSH configuration
            user_home_dirname (:obj:`str`, optional): user's home directory
        """
        ssh_config_filename = expand_user(
            os.path.join("~", ".ssh", "config"), user_home_dirname
        )

        with open(ssh_config_filename, "w") as file:
            file.write("".join(other_config).rstrip())

            file.write("\n\n")

            file.write("######## BEGIN DEEP ORIGIN CLI ########\n")
            file.write("".join(do_config).strip() + "\n")
            file.write("######## END DEEP ORIGIN CLI ########\n")
