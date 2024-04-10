from .env_var import EnvironmentVariable


class SecretEnvironmentVariable(EnvironmentVariable):
    """Secret environment variable"""

    class Meta:
        platform_id = "Secret environment variables"
