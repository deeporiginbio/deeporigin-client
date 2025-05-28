"""class to handle Anthropic API key"""

import pydantic

from .secret_env_var_value import SecretEnvironmentVariableValue


class AnthropicApiKey(SecretEnvironmentVariableValue):
    """Anthropic API key"""

    class Meta:
        platform_id = "Anthropic API key"

    @classmethod
    @property
    def KEY(cls) -> str:
        return "ANTHROPIC_API_KEY"

    value: str = pydantic.Field(default=None, min_length=1, validate_default=True)
