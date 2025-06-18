"""class to handle Anthropic API key"""

import pydantic

from .secret_env_var_value import SecretEnvironmentVariableValue


class AnthropicApiKey(SecretEnvironmentVariableValue):
    """Anthropic API key"""

    class Meta:
        """Meta class for AnthropicApiKey"""

        platform_id = "Anthropic API key"

    @classmethod
    def KEY(cls) -> str:
        """Anthropic API key"""

        return "ANTHROPIC_API_KEY"

    value: str = pydantic.Field(default=None, min_length=1, validate_default=True)
