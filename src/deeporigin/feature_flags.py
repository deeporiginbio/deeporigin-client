"""module to interact with feature flags"""

import dataclasses
import functools

from deeporigin.config import get_value as get_config
from deeporigin.warnings import DeepOriginWarning

__all__ = [
    "get_value",
    "FeatureNotAvailableWarning",
    "FeatureFlags",
]


@dataclasses.dataclass
class FeatureFlags:
    """Feature flags"""

    variables: bool


@functools.cache
def get_value() -> FeatureFlags:
    """Get the values of each feature flag

    Returns:
        :obj:`FeatureFlags`: values of the feature flags
    """
    config = get_config()

    if config.feature_flags:
        return FeatureFlags(**config.feature_flags)

    flags = FeatureFlags(variables=True)
    return flags


class FeatureNotAvailableWarning(DeepOriginWarning):
    """Warning raised when a user tries to access a feature that is not available to them"""

    pass
