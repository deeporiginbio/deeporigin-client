def pytest_addoption(parser):
    parser.addoption("--client", action="store", default="mock")
    parser.addoption("--responses", action="store", default="pass", help="Set to store to stash responses")


def pytest_generate_tests(metafunc):
    # Parametrize "client" fixture if it exists in the test function arguments
    option_value = metafunc.config.option.client
    if "client" in metafunc.fixturenames and option_value is not None:
        metafunc.parametrize("client", [option_value])

    option_value = metafunc.config.option.responses
    if "responses" in metafunc.fixturenames and option_value is not None:
        metafunc.parametrize("responses", [option_value])