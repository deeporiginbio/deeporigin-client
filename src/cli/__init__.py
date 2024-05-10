import functools
import sys
import warnings

import cement
import termcolor
from deeporigin import auth
from deeporigin.config import get_value
from deeporigin.utils import _print_dict

from .. import __version__
from ..exceptions import DeepOriginException
from ..warnings import DeepOriginWarning
from .context import CONTROLLERS as CONTEXT_CONTROLLERS
from .data import CONTROLLERS as DATA_CONTROLLERS
from .variables import CONTROLLERS as VARIABLE_CONTROLLERS

__all__ = ["main", "App"]


class BaseController(cement.Controller):
    class Meta:
        label = "base"
        help = "Utility for managing Deep Origin ComputeBenches"
        description = (
            "Utility for managing Deep Origin ComputeBenches, such as retrieving variables and "
            "secrets from the Deep Origin platform and installing them into benches."
        )
        arguments = [
            (
                ["-v", "--version"],
                {
                    "action": "store_true",
                    "help": "Display the version of this application",
                },
            ),
        ]

    @cement.ex(hide=True)
    def _default(self):
        args = self.app.pargs
        if args.version:
            print(__version__)
        else:
            raise SystemExit(self._parser.print_help())

    @cement.ex(
        help="Authenticate to the Deep Origin platform",
    )
    def authenticate(self):
        """list the columns of the row and their values, where applicable"""
        auth.get_tokens(refresh=False)

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
    def config(self):
        """list the columns of the row and their values, where applicable"""

        data = get_value()

        data.pop("list_bench_variables_query_template", None)

        _print_dict(data, json=self.app.pargs.json)


class App(cement.App):
    class Meta:
        label = "deeporigin"
        base_controller = "base"
        handlers = (
            [BaseController]
            + VARIABLE_CONTROLLERS
            + CONTEXT_CONTROLLERS
            + DATA_CONTROLLERS
        )


def except_hook(built_in_except_hook, type, value, tb):
    if issubclass(type, DeepOriginException):
        sys.stderr.write(termcolor.colored(value, "red") + "\n")
    else:
        built_in_except_hook(type, value, tb)


def set_highlighted_except_hook():
    built_in_except_hook = sys.excepthook

    sys.excepthook = functools.partial(except_hook, built_in_except_hook)


def format_warning(built_in_format_warning, msg, category, *args, **kwargs):
    if issubclass(category, DeepOriginWarning):
        return termcolor.colored(str(msg), "yellow") + "\n"
    else:
        return built_in_format_warning(msg, category, *args, **kwargs)


def set_format_warning():
    built_in_format_warning = warnings.formatwarning
    warnings.formatwarning = functools.partial(format_warning, built_in_format_warning)


def main():
    set_highlighted_except_hook()
    set_format_warning()
    with App() as app:
        app.run()


if __name__ == "__main__":
    main()
