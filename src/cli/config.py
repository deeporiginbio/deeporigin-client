"""This implements tools to configure the CLI and Python client"""

import os
import shutil

import cement
import yaml
from deeporigin.config import (
    CONFIG_YML_LOCATION,
    TEMPLATE,
    get_value,
)
from deeporigin.exceptions import DeepOriginException
from deeporigin.utils import _print_dict


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

        data.pop("list_workstation_variables_query_template", None)

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

        # check that key exists in the confuse template
        if key not in TEMPLATE.keys():
            raise DeepOriginException(
                message=f"{key} is not a valid key for a configuration file.",
                fix=f" Should be one of {', '.join(list(TEMPLATE.keys()))}",
            )

        # check if config file exists
        if os.path.isfile(CONFIG_YML_LOCATION):
            with open(CONFIG_YML_LOCATION, "r") as file:
                data = yaml.safe_load(file)
        else:
            # no file.
            data = {}
        data[key] = value

        # check that this is valid
        get_value(override_values=tuple(data.items()))

        # save configuration
        with open(CONFIG_YML_LOCATION, "w") as file:
            yaml.dump(data, file, default_flow_style=False)

        print(f"✔︎ {key} → {value}")

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
        if not os.path.isfile(CONFIG_YML_LOCATION):
            raise DeepOriginException(
                "Cannot save configuration for later use because no user configuration exists."
            )

        save_location = os.path.expanduser(f"~/.deeporigin/{name}.yml")
        shutil.copy(CONFIG_YML_LOCATION, save_location)
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

        file_to_load = os.path.expanduser(f"~/.deeporigin/{name}.yml")

        # check if config file exists
        if not os.path.isfile(file_to_load):
            raise DeepOriginException(
                "Cannot save configuration for later use because no user configuration exists."
            )

        shutil.copy(file_to_load, CONFIG_YML_LOCATION)
        print(f"✔︎ Configuration loaded from {file_to_load}")

        # loading a config should remove api_tokens to trigger a re-authentication
        try:
            os.remove(os.path.expanduser("~/.deeporigin/api_tokens"))
        except FileNotFoundError:
            pass


CONTROLLERS = [
    ConfigController,
]
