import os

from ... import utils
from ..base_type import Variable
from .env_var import EnvironmentVariable
from .secret_file_value import SecretFileValue


class XpressLicenseFile(SecretFileValue):
    """XPRESS license file"""

    class Meta:
        platform_id = "XPRESS license file"

    @classmethod
    @property
    def FILENAME(cls) -> str:
        return "~/xpauth.xpr"

    @classmethod
    @property
    def KEY(cls) -> str:
        return "XPAUTH_PATH"

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
        super().install(
            prev_self,
            user_home_dirname=user_home_dirname,
            overwrite=overwrite,
        )

        filename = utils.expand_user(self.FILENAME, user_home_dirname=user_home_dirname)
        EnvironmentVariable.install_env_var(
            self.name,
            self.KEY,
            filename,
            filename if prev_self else None,
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
        uninstalled = super().uninstall(
            user_home_dirname=user_home_dirname, overwrite=overwrite
        )

        filename = utils.expand_user(self.FILENAME, user_home_dirname=user_home_dirname)
        EnvironmentVariable.uninstall_env_var(
            self.name,
            self.KEY,
            filename,
            is_user_variable=False,
            user_home_dirname=user_home_dirname,
            overwrite=overwrite,
        )

        return uninstalled
