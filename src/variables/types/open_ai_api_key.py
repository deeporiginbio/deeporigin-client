import pydantic

from .secret_env_var_value import SecretEnvironmentVariableValue


class OpenAiApiKey(SecretEnvironmentVariableValue):
    """Open AI API key."""

    class Meta:
        platform_id = "Open AI API key"

    @classmethod
    @property
    def KEY(cls) -> str:
        return "OPENAI_API_KEY"

    value: str = pydantic.Field(default=None, min_length=1, validate_default=True)
