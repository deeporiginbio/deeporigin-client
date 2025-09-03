"""this module tests the tools API"""

import pytest

from deeporigin.platform import Client, tools_api
from tests.utils import config  # noqa: F401


def test_tools_api(config):  # noqa: F811
    """test the tools API"""

    if config["mock"]:
        pytest.skip("test skipped with mock client")

    # Create client inside the test function to avoid hanging at import time
    try:
        client = Client()
    except Exception as e:
        pytest.skip(f"Could not create client: {e}")

    tools = tools_api.get_all_tools(client=client)
    assert len(tools) > 0

    print(f"Found {len(tools)} tools")
