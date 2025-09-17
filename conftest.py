"""pytest configuration file

this governs the arguments we can pass pytest.

If you pass --mock, the tests will use a mock client.
If you pass --record, the tests will record responses (only useful for a live client)
"""


def pytest_addoption(parser):
    parser.addoption(
        "--mock", action="store_true", default=False, help="Use mock client"
    )
    parser.addoption(
        "--record", action="store_true", default=False, help="Set to record responses"
    )
    parser.addoption(
        "--org_key",
        action="store",
        default="deeporigin",
        help="Organization key to use for Client/MockClient",
    )


def pytest_generate_tests(metafunc):
    # Parametrize "mock" fixture if it exists in the test function arguments
    option_value = metafunc.config.option.mock
    if "mock" in metafunc.fixturenames and option_value is not None:
        metafunc.parametrize("mock", [option_value])

    option_value = metafunc.config.option.record
    if "record" in metafunc.fixturenames and option_value is not None:
        metafunc.parametrize("record", [option_value])

    option_value = getattr(metafunc.config.option, "org_key", None)
    if "org_key" in metafunc.fixturenames and option_value is not None:
        metafunc.parametrize("org_key", [option_value])
