"""this implements controllers and hooks to connect to
managed_data.py"""

import os

import cement
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
                fix=f" Should be one of {", ".join(list(TEMPLATE.keys()))}",
            )

        # check if config file exists
        if os.path.isfile(CONFIG_YML_LOCATION):
            print("need to read file")
            dsfsdf
        else:
            # no file.
            with open(CONFIG_YML_LOCATION, "w") as file:
                yaml.dump({key: value}, file, default_flow_style=False)


CONTROLLERS = [
    ConfigController,
]
