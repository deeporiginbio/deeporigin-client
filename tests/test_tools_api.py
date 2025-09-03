"""this module tests the tools API"""

import pytest

from deeporigin.platform import Client, tools_api
from tests.utils import config  # noqa: F401

try:
    client = Client()
except Exception:
    client = None


def test_tools_api(config):  # noqa: F811
    """test the tools API"""

    if config["mock"]:
        pytest.skip("test skipped with mock client")

    tools = tools_api.get_all_tools(client=client)
    assert len(tools) > 0

    print(f"Found {len(tools)} tools")
