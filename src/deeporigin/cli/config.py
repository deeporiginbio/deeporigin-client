"""This implements tools to configure the CLI and Python client"""

import os
import shutil

import cement
import yaml

from deeporigin.config import CONFIG_YML_LOCATION, get_value, set_value
from deeporigin.exceptions import DeepOriginException
from deeporigin.utils.core import (
    _ensure_do_folder,
    _get_api_tokens_filepath,
    _print_dict,
)


class ConfigController(cement.Controller):
    """Controller for the config subcommand of the CLI"""

    class Meta:
        label = "config"
        stacked_on = "base"
        stacked_type = "nested"
        help = "Show and modify configuration"
        description = """Show and modify the configuration for this application"""

    def _default(self):
        self.app.pargs.json = False
        self.show()

    @cement.ex(
        help="Show configuration",
        arguments=[
            (
                ["--json"],
                {
                    "action": "store_true",
                    "help": "Whether to return JSON formatted data [default: False]",
                },
            ),
        ],
    )
    def show(self):
        """Show configuration"""

        data = get_value()

        _print_dict(data, json=self.app.pargs.json, key_label="Variable")

    @cement.ex(
        help="Set configuration value",
        arguments=[
            (["key"], {"help": "Key to set", "action": "store"}),
            (["value"], {"help": "Value to set", "action": "store"}),
        ],
    )
    def set(self):
        """Set configuration value and save to config file"""

        key = self.app.pargs.key
        value = self.app.pargs.value

        set_value(key, value)

    @cement.ex(
        help="Save configuration to a file",
        arguments=[
            (["profile"], {"help": "Profile to save as", "action": "store"}),
        ],
    )
    def save(self):
        """Save configuration to a file"""

        name = self.app.pargs.profile

        # check if config file exists
        save_location = _ensure_do_folder() / f"{name}.yml"

        if os.path.isfile(CONFIG_YML_LOCATION):
            shutil.copy(CONFIG_YML_LOCATION, save_location)
        else:
            with open(save_location, "w") as file:
                yaml.dump({}, file, default_flow_style=False)

        print(f"✔︎ Configuration saved to {save_location}")

    @cement.ex(
        help="Load configuration from a file",
        arguments=[
            (["profile"], {"help": "Profile name to save as", "action": "store"}),
        ],
    )
    def load(self):
        """Load configuration from a file"""

        name = self.app.pargs.profile

        file_to_load = _ensure_do_folder() / f"{name}.yml"

        # check if config file exists
        if not os.path.isfile(file_to_load):
            raise DeepOriginException(
                f"Configuration '{name}' could not be loaded because it does not exist."
            )

        shutil.copy(file_to_load, CONFIG_YML_LOCATION)
        print(f"✔︎ Configuration loaded from {file_to_load}")

        # loading a config should remove api_tokens to trigger a re-authentication
        try:
            os.remove(_get_api_tokens_filepath())
        except FileNotFoundError:
            pass


CONTROLLERS = [
    ConfigController,
]
