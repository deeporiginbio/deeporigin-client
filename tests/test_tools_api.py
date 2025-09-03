"""this module tests the tools API"""

from deeporigin.platform import Client, tools_api

try:
    client = Client()
except Exception:
    client = None


def test_tools_api():
    """test the tools API"""
    tools = tools_api.get_all_tools(client=client)
    assert len(tools) > 0

    print(f"Found {len(tools)} tools")
