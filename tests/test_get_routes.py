"""this module tests various GET routes on the platform API. it is meant to be run against a live instance"""

import pytest

import deeporigin.platform.tools_api as tools_api
from tests.utils import config  # noqa: F401


def test_get_tool_executions(config):  # noqa: F811
    if config["mock"]:
        pytest.skip("test skipped with mock client")
    jobs = tools_api.get_tool_executions(
        org_key="deeporigin",
        filter={},
    )
    print(jobs)


def test_get_all_tools(config):  # noqa: F811
    if config["mock"]:
        pytest.skip("test skipped with mock client")
    tools = tools_api.get_all_tools()

    assert isinstance(tools, list)
    assert len(tools) > 0


def test_get_functions(config):  # noqa: F811
    if config["mock"]:
        pytest.skip("test skipped with mock client")
    functions = tools_api.get_functions()

    assert isinstance(functions, list)
    assert len(functions) > 0


def test_get_executions(config):  # noqa: F811
    if config["mock"]:
        pytest.skip("test skipped with mock client")
    jobs = tools_api.get_tool_executions()
    assert isinstance(jobs, list)
    assert len(jobs) > 0


def test_get_functions_bykey(config):  # noqa: F811
    if config["mock"]:
        pytest.skip("test skipped with mock client")
    functions = tools_api.get_functions()

    assert isinstance(functions, list)
    assert len(functions) > 0

    functions = tools_api.get_functions_bykey(key=functions[0].key)

    assert isinstance(functions, list)


def test_get_job_df(config):  # noqa: F811
    if config["mock"]:
        pytest.skip("test skipped with mock client")
    import pandas as pd

    from deeporigin.tools.job import get_dataframe

    df = get_dataframe()
    assert isinstance(df, pd.DataFrame)
