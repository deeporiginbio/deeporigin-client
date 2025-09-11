"""this module tests the tools API"""

import time

import pytest

from deeporigin.platform import Client, tools_api
from tests.utils import config  # noqa: F401


@pytest.mark.dependency()
def test_health(config):  # noqa: F811
    """test the health API"""

    if config["mock"]:
        pytest.skip("test skipped with mock client")

    data = tools_api.check(client=Client())
    assert data["status"] == "ok"
    assert data["info"]["mikroOrm"]["status"] == "up"


@pytest.mark.dependency(depends=["test_health"])
def test_get_all_tools(config):  # noqa: F811
    """test the tools API"""

    if config["mock"]:
        pytest.skip("test skipped with mock client")

    tools = tools_api.get_all_tools(client=Client())
    assert len(tools) > 0

    print(f"Found {len(tools)} tools")


@pytest.mark.dependency(depends=["test_health"])
def test_get_all_functions(config):  # noqa: F811
    """test the functions API"""

    if config["mock"]:
        pytest.skip("test skipped with mock client")

    functions = tools_api.get_functions(client=Client())
    assert len(functions) > 0

    print(f"Found {len(functions)} functions")


@pytest.mark.dependency(depends=["test_health"])
def test_get_all_executions(config):  # noqa: F811
    """test the executions API"""

    if config["mock"]:
        pytest.skip("test skipped with mock client")

    executions = tools_api.get_tool_executions(client=Client())

    print(f"Found {len(executions)} executions")


def test_job(config):  # noqa: F811
    if config["mock"]:
        pytest.skip("test skipped with mock client")

    from deeporigin.tools.job import Job

    jobs = tools_api.get_tool_executions()
    execution_id = jobs[0].executionId
    job = Job.from_id(execution_id)

    assert execution_id in job._ids


def test_job_df(config):  # noqa: F811
    if config["mock"]:
        pytest.skip("test skipped with mock client")

    from deeporigin.tools.job import get_dataframe

    _ = get_dataframe()


@pytest.mark.dependency()
def test_job_df_filtering(config):  # noqa: F811
    if config["mock"]:
        pytest.skip("test skipped with mock client")

    from deeporigin.drug_discovery.constants import tool_mapper
    from deeporigin.tools.job import get_dataframe

    tool_key = tool_mapper["Docking"]

    df = get_dataframe(
        tool_key=tool_key,
    )

    assert len(df["tool_key"].unique()) == 1, (
        f"should only be one tool key. Instead there were {len(df['tool_key'].unique())}"
    )

    assert df["tool_key"].unique()[0] == tool_key, (
        f"Expected to get back jobs for {tool_key}. Instead got {df['tool_key'].unique()[0]}"
    )


@pytest.mark.dependency(depends=["test_job_df_filtering"])
def test_run_docking_and_cancel(config):  # noqa: F811
    if config["mock"]:
        pytest.skip("test skipped with mock client")

    from deeporigin.drug_discovery import (
        BRD_DATA_DIR,
        Complex,
    )

    sim = Complex.from_dir(BRD_DATA_DIR)

    pockets = sim.protein.find_pockets(pocket_count=1, use_cache=False)
    pocket = pockets[0]

    job = sim.docking.run(
        pocket=pocket,
        re_run=True,
        use_parallel=False,
        client=Client(),
    )

    # wait for a bit to start
    time.sleep(10)

    # check that it's running
    job.sync()
    assert "Running" in job._status

    # now cancel it
    job.cancel()

    # wait for a bit to cancel
    time.sleep(10)

    # check that it's cancelled
    job.sync()
    assert "Cancelled" in job._status
