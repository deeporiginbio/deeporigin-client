import cement

from ..context import (
    get_value as get_context,
)

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
        help = "Get the context for the Deep Origin ComputeBench"
        description = "Get the context for the Deep Origin ComputeBench, such as the ID of the bench, user and organization, host compute cluster, and hardware blueprint."
        arguments = []

    @cement.ex(hide=True)
    def _default(self):
        """Default action. returns the context"""
        context = get_context()
        print(context)


CONTROLLERS = [
    ContextController,
]
