"""this module tests the tools API"""

import pytest

from deeporigin.platform import Client, tools_api
from tests.utils import config  # noqa: F401


def test_get_all_tools(config):  # noqa: F811
    """test the tools API"""

    if config["mock"]:
        pytest.skip("test skipped with mock client")

    tools = tools_api.get_all_tools(client=Client())
    assert len(tools) > 0

    print(f"Found {len(tools)} tools")


def test_get_all_functions(config):  # noqa: F811
    """test the functions API"""

    if config["mock"]:
        pytest.skip("test skipped with mock client")

    functions = tools_api.get_functions(client=Client())
    assert len(functions) > 0

    print(f"Found {len(functions)} functions")


def test_get_all_executions(config):  # noqa: F811
    """test the executions API"""

    if config["mock"]:
        pytest.skip("test skipped with mock client")

    executions = tools_api.get_tool_executions(client=Client())

    print(f"Found {len(executions)} executions")
