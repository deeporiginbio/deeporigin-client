"""class to handle OpenAI API keys"""

import pydantic

from .secret_env_var_value import SecretEnvironmentVariableValue


class OpenAiApiKey(SecretEnvironmentVariableValue):
    """Open AI API key."""

    class Meta:
        """Meta class for OpenAiApiKey"""

        platform_id = "Open AI API key"

    @classmethod
    def KEY(cls) -> str:
        """Open AI API key"""

        return "OPENAI_API_KEY"

    value: str = pydantic.Field(default=None, min_length=1, validate_default=True)
