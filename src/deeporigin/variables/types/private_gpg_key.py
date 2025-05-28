"""class to handle GPG keys"""

import os
import subprocess
import tempfile
import typing

import pydantic

from ...exceptions import DeepOriginException
from ..base_type import Variable

T = typing.TypeVar("T")


class PrivateGpgKey(Variable):
    """Private GPG key"""

    class Meta:
        platform_id = "GPG signing keys"

    value: str = pydantic.Field(validate_default=True)

    @pydantic.field_validator("value")
    @classmethod
    def value_must_be_a_valid_private_gpg_key(cls, value: str) -> str:
        """Verify whether a value is a valid GPG private key

        Args:
            value (:obj:`str`): value

        Returns:
            :obj:`str`: value

        Raises:
            :obj:`ValueError`: if the value is not a valid GPG private key
        """
        if not isinstance(value, str):
            raise ValueError(
                f"Value must be a GPG private key. {value} is not a GPG private key."
            )

        if not (value.startswith("-----BEGIN PGP PRIVATE KEY BLOCK-----")):
            raise ValueError(
                f"Value must be a GPG private key. {value} is not a GPG private key."
            )

        _, key_filename = tempfile.mkstemp()
        with open(key_filename, "w") as file:
            file.write(value)

        user_home_dirname = os.path.expanduser("~")
        gpg_dirname = os.path.join(user_home_dirname, ".gnupg")

        cmd = []
        cmd.extend(["mkdir", "-p", gpg_dirname])

        cmd.append("&&")
        cmd.extend(["chmod", "700", gpg_dirname])
        cmd.append("&&")

        cmd.extend(["gpg", "--batch", "--import", "--dry-run", key_filename])

        completed_process = subprocess.run(
            " ".join(cmd),
            shell=True,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )

        os.remove(key_filename)

        if completed_process.returncode not in [0, 130]:
            raise ValueError(
                f"Value must be a GPG private key. {value} is not a GPG private key."
            )

        return value

    def install(
        self,
        prev_self: Variable,
        user_home_dirname: str = os.path.expanduser("~"),
        overwrite: bool = False,
    ) -> None:
        """Install a variable/secret: Run `gpg --import` with the key

        Args:
            prev_self (:obj:`Variable`): previous value of the variable
            user_home_dirname (:obj:`str`, optional): user's home directory
            overwrite (:obj:`bool`, optional): Whether to overwrite direct changes to variables and secrets
        """
        value = self.value

        _, key_filename = tempfile.mkstemp()
        with open(key_filename, "w") as file:
            file.write(value)

        cmd = ["gpg", "--batch", "--import", key_filename]
        completed_process = subprocess.run(
            cmd,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )

        os.remove(key_filename)

        if completed_process.returncode not in [0, 130]:
            raise DeepOriginException(
                f"Private GPG key {self.name or ''} could not be imported: {completed_process.stdout.decode().strip()}"
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
        value = self.value

        _, key_filename = tempfile.mkstemp()
        with open(key_filename, "w") as file:
            file.write(value)

        cmd = ["gpg", "--show-keys", key_filename]
        completed_process = subprocess.run(
            cmd,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )

        os.remove(key_filename)

        stdout = completed_process.stdout.decode()

        if completed_process.returncode not in [0, 130]:
            raise DeepOriginException(
                f"Private GPG key {self.name or ''} could not be inspected: {stdout.strip()}"
            )

        uninstalled = True

        key_id = stdout.partition("\n")[2].partition("\n")[0].strip()

        cmd = [
            "expect",
            "-c",
            f'spawn gpg --batch --delete-secret-keys {key_id}; expect "Do you really want" {{send "\\r"}}; interact',
        ]
        completed_process = subprocess.run(
            cmd,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )

        if completed_process.returncode not in [0, 130]:
            stdout = completed_process.stdout.decode().strip()
            raise DeepOriginException(
                f"Private GPG key {self.name or ''} could not be uninstalled: {stdout}"
            )

        cmd = ["gpg", "--batch", "--delete-keys", "--fingerprint", key_id]
        completed_process = subprocess.run(
            cmd,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )

        if completed_process.returncode not in [0, 130]:
            stdout = completed_process.stdout.decode().strip()
            raise DeepOriginException(
                f"Private GPG key {self.name or ''} could not be uninstalled: {stdout}"
            )

        return uninstalled
