def pytest_addoption(parser):
    parser.addoption("--client", action="store", default="mock")


def pytest_generate_tests(metafunc):
    # This is called for every test. Only get/set command line arguments
    # if the argument is specified in the list of test "fixturenames".
    option_value = metafunc.config.option.client
    if "client" in metafunc.fixturenames and option_value is not None:
        metafunc.parametrize("client", [option_value])
