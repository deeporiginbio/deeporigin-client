"""helper module to set up tests"""

from mock_client import MockClient
import pytest

from deeporigin.platform import Client


@pytest.fixture(scope="session", autouse=True)
def config(pytestconfig):
    """this fixture performs some setup tasks
    before all tests are run, and runs only once"""

    data = {}

    # set up client
    if pytestconfig.getoption("mock"):
        client = MockClient()
        data["mock"] = True

    else:
        client = Client()
        data["mock"] = False
        if pytestconfig.getoption("record"):
            client.recording = True

    data["client"] = client

    # tests run on yield
    yield data
