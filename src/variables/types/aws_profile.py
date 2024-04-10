import os
import re
import subprocess
import typing
import warnings

import pydantic

from ...exceptions import DeepOriginException
from ...warnings import DeepOriginWarning
from ..base_type import Variable

T = typing.TypeVar("T")


class AwsProfile(Variable):
    """AWS profile"""

    class Meta:
        platform_id = "AWS Credentials"

    profile: str = pydantic.Field(min_length=1, validate_default=True)
    access_key_id: str = pydantic.Field(
        pattern=r"^[A-Z0-9]{20}$", validate_default=True
    )
    secret_access_key: str = pydantic.Field(
        pattern=r"^[A-Za-z0-9/+=]{40}$", validate_default=True
    )
    region: typing.Optional[str] = pydantic.Field(
        pattern=r"^[a-z]+\-[a-z]+\-\d+[a-z]?$", default=None, validate_default=True
    )

    @classmethod
    def from_platform(cls: type[T], platform_value: dict) -> T:
        """Create an instance from a dictionary from the Deep Origin platform

        Args:
            platform_value (:obj:`dict`): platform representation of the variable/secret

        Returns:
            :obj:`AwsProfile`: Python instance of the variable/secret

        Raises:
            :obj:`DeepOriginException`: if value is invalid
        """
        return cls(
            drn=platform_value["drn"],
            name=platform_value["name"],
            profile=platform_value["profileName"],
            access_key_id=platform_value["accessKey"],
            secret_access_key=platform_value["secretKey"],
            region=platform_value["defaultRegion"] or None,
        )

    def install(
        self,
        prev_self: Variable,
        user_home_dirname: str = os.path.expanduser("~"),
        overwrite: bool = False,
    ) -> None:
        """Install a variable/secret: Append profile to '~/.aws/credentials`

        Args:
            prev_self (:obj:`Variable`): previous value of the variable
            user_home_dirname (:obj:`str`, optional): user's home directory
            overwrite (:obj:`bool`, optional): Whether to overwrite direct changes to variables and secrets
        """
        if overwrite or prev_self is None:
            prev_access_key_id = None
            prev_secret_access_key = None
            prev_region = None
        else:
            (
                prev_access_key_id,
                prev_secret_access_key,
                prev_region,
            ) = self.get_aws_profile(self.profile)

        if (
            overwrite
            or prev_self is None
            or (
                prev_access_key_id == prev_self.access_key_id
                and prev_secret_access_key == prev_self.secret_access_key
                and prev_region == prev_self.region
            )
        ):
            access_key_id = self.access_key_id
            secret_access_key = self.secret_access_key
            region = self.region

            cmd = []

            if access_key_id:
                cmd.append("&&")
                cmd.extend(
                    [
                        "aws",
                        "configure",
                        "--profile",
                        self.profile,
                        "set",
                        "aws_access_key_id",
                        access_key_id,
                    ]
                )

            if secret_access_key:
                cmd.append("&&")
                cmd.extend(
                    [
                        "aws",
                        "configure",
                        "--profile",
                        self.profile,
                        "set",
                        "aws_secret_access_key",
                        secret_access_key,
                    ]
                )

            if region:
                cmd.append("&&")
                cmd.extend(
                    [
                        "aws",
                        "configure",
                        "--profile",
                        self.profile,
                        "set",
                        "region",
                        region,
                    ]
                )

            if cmd:
                completed_process = subprocess.run(
                    " ".join(cmd[1:]),
                    shell=True,
                    check=False,
                    stdout=subprocess.PIPE,
                    stderr=subprocess.STDOUT,
                )

                if completed_process.returncode != 0:
                    raise DeepOriginException(
                        f"AWS profile {self.name or ''} could not be imported: {completed_process.stdout.decode().strip()}"
                    )
        else:
            warnings.warn(
                (
                    f"Variable {self.name or ''} was not modified because it has been modified directly. "
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
        prev_access_key_id, prev_secret_access_key, prev_region = self.get_aws_profile(
            self.profile
        )

        uninstalled = False

        if overwrite or (
            prev_access_key_id == self.access_key_id
            and prev_secret_access_key == self.secret_access_key
            and prev_region == self.region
        ):
            uninstalled = True

            # TODO: improve parsing of AWS credentials and profiles
            cmd = []

            aws_credentials_filename = os.path.join(
                user_home_dirname, ".aws", "credentials"
            )
            if os.path.isfile(aws_credentials_filename):
                sub = re.escape(self.profile).replace("'", "'")
                cmd.extend(
                    [
                        "sed",
                        "--in-place",
                        "--regexp-extended",
                        rf"'/\[{sub}\]/d'",
                        aws_credentials_filename,
                    ]
                )

                cmd.append("&&")
                sub = re.escape(self.access_key_id).replace("'", "'")
                cmd.extend(
                    [
                        "sed",
                        "--in-place",
                        "--regexp-extended",
                        f"'/aws_access_key_id = {sub}/d'",
                        aws_credentials_filename,
                    ]
                )

                cmd.append("&&")
                sub = re.escape(self.secret_access_key).replace("'", "'")
                cmd.extend(
                    [
                        "sed",
                        "--in-place",
                        "--regexp-extended",
                        f"'/aws_secret_access_key = {sub}/d'",
                        aws_credentials_filename,
                    ]
                )

            # aws_config_filename = os.path.join(user_home_dirname, ".aws", "config")
            # if os.path.isfile(aws_config_filename):
            #     cmd.append("&&")
            #     sub = re.escape(self.profile).replace("'", "\'")
            #     cmd.extend(
            #         [
            #             "sed",
            #             "--in-place",
            #             "--regexp-extended",
            #             f"'/\[profile {sub}\]/d'",
            #             aws_config_filename,
            #         ]
            #     )
            #
            #     cmd.append("&&")
            #     sub = re.escape(self.region).replace("'", "\'")
            #     cmd.extend(
            #         [
            #             "sed",
            #             "--in-place",
            #             "--regexp-extended",
            #             f"'/region = {sub}/d'",
            #             aws_config_filename,
            #        ]
            #     )

            cmds = " ".join(cmd)
            completed_process = subprocess.run(
                cmds,
                shell=True,
                check=False,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )

            if completed_process.returncode != 0:
                raise DeepOriginException(
                    f"AWS profile {self.name or ''} could not be uninstalled: {completed_process.stdout.decode().strip()}"
                )
        else:
            warnings.warn(
                (
                    f"Variable {self.name or ''} was not uninstalled because it has been modified directly. "
                    "To update its value, rerun with the `overwrite` option."
                ),
                DeepOriginWarning,
            )

        return uninstalled

    @classmethod
    def get_aws_profile(cls, profile: str) -> tuple[str]:
        """Get the access key, secret key, and default region for an AWS

        Args:
            profile (:obj:`str`): profile name

        Returns:
            :obj:`tuple`: access key, secret key, and default region for the AWS profile
        """
        cmd = [
            "aws",
            "configure",
            "--profile",
            profile,
            "get",
            "aws_access_key_id",
        ]
        completed_process = subprocess.run(
            cmd,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        if completed_process.returncode == 0:
            access_key_id = completed_process.stdout.decode().strip()
        elif completed_process.returncode == 255:
            access_key_id = None
        else:
            raise DeepOriginException(
                f"AWS profile {profile} could not be inspected: {completed_process.stdout.decode().strip()}"
            )

        cmd = [
            "aws",
            "configure",
            "--profile",
            profile,
            "get",
            "aws_secret_access_key",
        ]
        completed_process = subprocess.run(
            cmd,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        if completed_process.returncode == 0:
            secret_access_key = completed_process.stdout.decode().strip()
        elif completed_process.returncode == 255:
            secret_access_key = None
        else:
            raise DeepOriginException(
                f"AWS profile {profile} could not be inspected: {completed_process.stdout.decode().strip()}"
            )

        cmd = [
            "aws",
            "configure",
            "--profile",
            profile,
            "get",
            "region",
        ]
        completed_process = subprocess.run(
            cmd,
            check=False,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
        )
        if completed_process.returncode == 0:
            region = completed_process.stdout.decode().strip()
        elif completed_process.returncode == 255:
            region = None
        else:
            raise DeepOriginException(
                f"AWS profile {profile} could not be inspected: {completed_process.stdout.decode().strip()}"
            )

        return access_key_id, secret_access_key, region
