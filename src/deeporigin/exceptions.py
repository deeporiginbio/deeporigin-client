"""custom exceptions to surface better errors"""

import sys
import textwrap
from typing import Optional

from tabulate import tabulate
import termcolor

__all__ = ["DeepOriginException"]


def _wrap(text: str) -> str:
    wrapped_message = textwrap.wrap(text)
    wrapped_message = "\n".join(wrapped_message)

    return wrapped_message


class DeepOriginException(Exception):
    """Deep Origin exception"""

    def __init__(
        self,
        message: str,
        *,
        title: str = "Deep Origin error",
        fix: Optional[str] = None,
    ):
        """Utility function to print a nicely formatted error, used in the CLI"""

        self.message = message

        printout = [[title], [_wrap(message)]]
        if fix:
            printout.append([_wrap(fix)])

        sys.stderr.write(
            termcolor.colored(
                tabulate(
                    printout,
                    tablefmt="rounded_grid",
                ),
                "red",
            )
            + "\n"
        )

        super().__init__(self.message)
