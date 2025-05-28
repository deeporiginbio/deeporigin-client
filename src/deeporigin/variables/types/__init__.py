import enum

from .anthropic_api_key import AnthropicApiKey
from .aws_profile import AwsProfile
from .env_var import EnvironmentVariable
from .file import File
from .git_http_credentials import GitHttpCredentials
from .gurobi_license_file import GurobiLicenseFile
from .mosek_license_file import MosekLicenseFile
from .open_ai_api_key import OpenAiApiKey
from .private_gpg_key import PrivateGpgKey
from .private_ssh_key import PrivateSshKey
from .secret_env_var import SecretEnvironmentVariable
from .secret_file import SecretFile
from .xpress_license_file import XpressLicenseFile

__all__ = [
    "EnvironmentVariable",
    "SecretEnvironmentVariable",
    "File",
    "SecretFile",
    "PrivateSshKey",
    "PrivateGpgKey",
    "GitHttpCredentials",
    "AwsProfile",
    "OpenAiApiKey",
    "AnthropicApiKey",
    "GurobiLicenseFile",
    "MosekLicenseFile",
    "XpressLicenseFile",
    "VariableType",
]


class VariableType(enum.Enum):
    """Type of variable/secret"""

    EnvironmentVariable = EnvironmentVariable
    SecretEnvironmentVariable = SecretEnvironmentVariable
    File = File
    SecretFile = SecretFile
    PrivateSshKey = PrivateSshKey
    PrivateGpgKey = PrivateGpgKey
    AwsProfile = AwsProfile
    GitHttpCredentials = GitHttpCredentials
    OpenAiApiKey = OpenAiApiKey
    AnthropicApiKey = AnthropicApiKey
    GurobiLicenseFile = GurobiLicenseFile
    MosekLicenseFile = MosekLicenseFile
    XpressLicenseFile = XpressLicenseFile
