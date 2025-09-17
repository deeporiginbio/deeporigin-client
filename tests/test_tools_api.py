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
    assert len(tools) > 0

    print(f"Found {len(tools)} tools")


@pytest.mark.dependency(depends=["test_tools_api_health"])
def test_get_all_functions(config):  # noqa: F811
    """test the functions API"""

    functions = tools_api.get_functions(client=config["client"])
    assert len(functions) > 0

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

    pockets = sim.protein.find_pockets(pocket_count=1, use_cache=True)
    pocket = pockets[0]

    job = sim.docking.run(
        pocket=pocket,
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
