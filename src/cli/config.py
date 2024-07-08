"""this implements controllers and hooks to connect to
managed_data.py"""

import os
import shutil

import cement
import confuse
import yaml
from deeporigin.config import CONFIG_YML_LOCATION, TEMPLATE, get_value
from deeporigin.exceptions import DeepOriginException
from deeporigin.utils import _print_dict


class ConfigController(cement.Controller):
    """Controller for data subcommand of CLI"""

    class Meta:
        label = "config"
        stacked_on = "base"
        stacked_type = "nested"
        help = "Show and modify configuration"
        description = """
Show and modify configuration file to connect to Deep Origin 
            """

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
        """list the columns of the row and their values, where applicable"""

        data = get_value()

        data.pop("list_bench_variables_query_template", None)

        _print_dict(data, json=self.app.pargs.json)

    @cement.ex(
        help="Set configuration value in config file",
        arguments=[
            (["key"], {"help": "Key to set", "action": "store"}),
            (["value"], {"help": "Value to set", "action": "store"}),
        ],
    )
    def set(self):
        """download or upload files or databases"""

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
            data[key] = value
        else:
            # no file.
            data = {key: value}

        # check that this is valid
        config = confuse.Configuration("deep_origin", __name__)
        config.set(data)
        try:
            _ = config.get(TEMPLATE)
        except confuse.confuse.ConfigError:
            raise DeepOriginException(
                f"{key}:{value} is not a valid entry for the configuration.",
            )

        with open(CONFIG_YML_LOCATION, "w") as file:
            yaml.dump(data, file, default_flow_style=False)

        print(f"✔︎ {key} → {value}")

    @cement.ex(
        help="Save configuration file for later use",
        arguments=[
            (["profile"], {"help": "Profile to save as", "action": "store"}),
        ],
    )
    def save(self):
        """save a config file for later use"""

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
        help="Save configuration file for later use",
        arguments=[
            (["profile"], {"help": "Profile name to save as", "action": "store"}),
        ],
    )
    def load(self):
        """load a config file and use it"""

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
