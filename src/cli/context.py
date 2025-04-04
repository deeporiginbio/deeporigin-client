"""controllers for context information in the CLI"""

import cement

from deeporigin.context import get_value as get_context

__all__ = [
    "CONTROLLERS",
    "ContextController",
]


class ContextController(cement.Controller):
    """Context controller"""

    class Meta:
        label = "context"
        stacked_on = "base"
        stacked_type = "nested"
        help = "Get information about a Deep Origin workstation"
        description = "Get information about a Deep Origin workstation, such as the ID of the workstation, user and organization, host compute cluster, and hardware."
        arguments = []

    @cement.ex(hide=True)
    def _default(self):
        """Default action. returns the context"""
        context = get_context()
        print(context)


CONTROLLERS = [
    ContextController,
]
