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
    org_key = pytestconfig.getoption("org_key")
    if pytestconfig.getoption("mock"):
        client = MockClient(org_key=org_key)
        data["mock"] = True

    else:
        client = Client(org_key=org_key) if org_key else Client()
        data["mock"] = False
        if pytestconfig.getoption("record"):
            client.recording = True

    data["client"] = client

    # tests run on yield
    yield data
