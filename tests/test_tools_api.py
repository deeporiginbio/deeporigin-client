"""this module tests the tools API"""

import time

import pytest

from deeporigin.drug_discovery import (
    BRD_DATA_DIR,
    Complex,
)
from deeporigin.drug_discovery.constants import tool_mapper
from deeporigin.platform import tools_api
from deeporigin.tools.job import Job, get_dataframe
from tests.utils import config  # noqa: F401

"""this module tests various GET routes on the platform API. it is meant to be run against a live instance"""


def test_get_tool_executions(config):  # noqa: F811
    jobs = tools_api.get_tool_executions(
        org_key="deeporigin",
        filter={},
        client=config["client"],
    )

    assert isinstance(jobs, list), "Expected a list"
    assert len(jobs) > 0, "Expected at least one job"


def test_get_functions(config):  # noqa: F811
    functions = tools_api.get_functions(client=config["client"])

    assert isinstance(functions, list), "Expected a list"
    assert len(functions) > 0, "Expected at least one function"


def test_get_executions(config):  # noqa: F811
    jobs = tools_api.get_tool_executions(client=config["client"])
    assert isinstance(jobs, list), "Expected a list"
    assert len(jobs) > 0, "Expected at least one job"


def test_get_functions_bykey(config):  # noqa: F811
    functions = tools_api.get_functions(client=config["client"])

    assert isinstance(functions, list), "Expected a list"
    assert len(functions) > 0, "Expected at least one function"

    functions = tools_api.get_functions_bykey(
        key=functions[0].key,
        client=config["client"],
    )

    assert isinstance(functions, list), "Expected a list"
    assert len(functions) > 0, "Expected at least one function"


def test_get_job_df(config):  # noqa: F811
    import pandas as pd

    df = get_dataframe(client=config["client"])
    assert isinstance(df, pd.DataFrame), "Expected a dataframe"
    assert len(df) > 0, "Expected at least one job"


@pytest.mark.dependency()
def test_tools_api_health(config):  # noqa: F811
    """test the health API"""

    data = tools_api.check(client=config["client"])
    assert data["status"] == "ok"
    assert data["info"]["mikroOrm"]["status"] == "up"


@pytest.mark.dependency(depends=["test_tools_api_health"])
def test_get_all_tools(config):  # noqa: F811
    """test the tools API"""

    tools = tools_api.get_all_tools(client=config["client"])
    assert len(tools) > 0, "Expected at least one tool"

    print(f"Found {len(tools)} tools")


@pytest.mark.dependency(depends=["test_tools_api_health"])
def test_get_all_functions(config):  # noqa: F811
    """test the functions API"""

    functions = tools_api.get_functions(client=config["client"])
    assert len(functions) > 0, "Expected at least one function"

    print(f"Found {len(functions)} functions")


@pytest.mark.dependency(depends=["test_tools_api_health"])
def test_get_all_executions(config):  # noqa: F811
    """test the executions API"""

    executions = tools_api.get_tool_executions(client=config["client"])

    print(f"Found {len(executions)} executions")


def test_job(config):  # noqa: F811
    jobs = tools_api.get_tool_executions(client=config["client"])
    execution_id = jobs[0].executionId
    job = Job.from_id(execution_id, client=config["client"])

    assert execution_id in job._ids


def test_job_df(config):  # noqa: F811
    _ = get_dataframe(client=config["client"])


@pytest.mark.dependency()
def test_job_df_filtering(config):  # noqa: F811
    tool_key = tool_mapper["Docking"]

    df = get_dataframe(
        tool_key=tool_key,
        client=config["client"],
    )

    assert len(df["tool_key"].unique()) == 1, (
        f"should only be one tool key. Instead there were {len(df['tool_key'].unique())}"
    )

    assert df["tool_key"].unique()[0] == tool_key, (
        f"Expected to get back jobs for {tool_key}. Instead got {df['tool_key'].unique()[0]}"
    )


def test_run_docking_and_cancel(config):  # noqa: F811
    sim = Complex.from_dir(BRD_DATA_DIR)
    sim.client = config["client"]

    job = sim.docking.run(
        box_size=(14.094597464129786, 14.094597464129786, 14.094597464129786),
        pocket_center=(-13.215283393859863, -6.083978652954102, 14.214159965515137),
        re_run=True,
        use_parallel=False,
    )

    # wait for a bit to start
    if not config["mock"]:
        time.sleep(10)

    # check that it's running
    job.sync()
    assert "Running" in job._status, f"Job with ID {job._ids} is not running"

    # now cancel it
    job.cancel()

    # wait for a bit to cancel
    if not config["mock"]:
        time.sleep(10)

    # check that it's cancelled
    job.sync()
    assert "Cancelled" in job._status, f"Job with ID {job._ids} is not cancelled"
