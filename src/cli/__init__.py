import functools
import sys
import warnings

import cement
import termcolor

from deeporigin import __version__, auth
from deeporigin.cli.config import CONTROLLERS as CONFIG_CONTROLLERS
from deeporigin.exceptions import DeepOriginException
from deeporigin.warnings import DeepOriginWarning

__all__ = ["main", "App"]


class BaseController(cement.Controller):
    class Meta:
        label = "base"
        help = "Client for Deep Origin"
        description = (
            "Client for Deep Origin such as for downloading data and installing variables and "
            "secrets into workstations."
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
        help="Authenticate to Deep Origin",
    )
    def authenticate(self):
        """list the columns of the row and their values, where applicable"""
        auth.get_tokens(refresh=False)


class App(cement.App):
    class Meta:
        label = "deeporigin"
        base_controller = "base"
        handlers = [BaseController] + CONFIG_CONTROLLERS


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
