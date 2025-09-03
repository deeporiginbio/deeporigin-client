"""Simplified configuration management for Deep Origin CLI/client.

This module stores and retrieves only two configuration values:
`env` and `org_key`.

Behavior:
- If the config file does not exist, it is created with `env=prod` and an
  empty `org_key`.
- If the config file exists, it is read and a dictionary is returned.
"""

import os

import yaml

from deeporigin.utils.core import _ensure_do_folder

CONFIG_YML_LOCATION = _ensure_do_folder() / "config.yml"

__all__ = ["get_value", "set_value", "CONFIG_YML_LOCATION"]


def _ensure_config_file_exists() -> None:
    """Ensure the configuration file exists; create with defaults if missing."""

    if not os.path.isfile(CONFIG_YML_LOCATION):
        default_data: dict = {"env": "prod", "org_key": ""}
        os.makedirs(os.path.dirname(CONFIG_YML_LOCATION), exist_ok=True)
        with open(CONFIG_YML_LOCATION, "w") as file:
            yaml.safe_dump(default_data, file, default_flow_style=False)


def get_value() -> dict:
    """Get the configuration values.

    Creates the file with defaults if it doesn't exist, then returns a dict
    with keys `env` and `org_key`.

    Args:
        config_file_location: Optional custom path for the config file.

    Returns:
        A dictionary with keys `env` and `org_key`.
    """

    _ensure_config_file_exists()

    with open(CONFIG_YML_LOCATION, "r") as file:
        data = yaml.safe_load(file) or {}

    # Fill defaults if missing
    env = data.get("env", "prod")
    org_key = data.get("org_key", "")

    return {"env": env, "org_key": org_key}


def set_value(key: str, value) -> None:
    """Set a configuration value.

    Only `env` and `org_key` are supported keys.

    Args:
        key: Configuration key to set (must be `env` or `org_key`).
        value: Value to set.
    """

    if key not in {"env", "org_key"}:
        raise ValueError(
            f"{key} is not a valid configuration key. Supported keys are: env, org_key"
        )

    _ensure_config_file_exists()

    with open(CONFIG_YML_LOCATION, "r") as file:
        data = yaml.safe_load(file) or {}

    data[key] = value

    # Persist updated data
    with open(CONFIG_YML_LOCATION, "w") as file:
        yaml.safe_dump(data, file, default_flow_style=False)

    print(f"✔︎ {key} → {value}")
