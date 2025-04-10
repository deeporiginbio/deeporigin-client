import json
from typing import Union

from beartype import beartype
import pytest

from tests.utils import (
    _run_cli_command,
    config,  # noqa: F401
)

# this allows us to try every CLI command with both
# multiple options using pytest.mark.parametrize
JSON_OPTIONS = (
    ["--json"],
    [],
)
LIST_OPTIONS = (
    [],
    ["--rows"],
    ["--databases"],
    ["--files"],
    ["--folders"],
    ["--folders", "--databases"],
    ["--folders", "--rows"],
    ["--databases", "--rows"],
)


def test_data(config):  # noqa: F811
    stdout = _run_cli_command(
        ["data"],
        client=config["client"],
    )

    assert "List data in the data hub on Deep Origin" in stdout, "Unexpected output"


@pytest.mark.parametrize("json_option", JSON_OPTIONS)
def test_describe_file(config, json_option):  # noqa: F811
    file_id = config["file"].id

    stdout = _run_cli_command(
        ["data", "describe", file_id] + json_option,
        client=config["client"],
    )

    if json_option == ["--json"]:
        data = _check_json(stdout)
        assert data["id"] == file_id, "Expected file ID to match"


@pytest.mark.parametrize("option", JSON_OPTIONS)
def test_describe_row(config, option):  # noqa: F811
    row_id = config["rows"][0]

    stdout = _run_cli_command(
        ["data", "describe", row_id] + option,
        client=config["client"],
    )
    assert row_id in stdout, "Expected to see row_id in output"

    # check that we can parse into JSON
    if option == ["--json"]:
        data = _check_json(stdout)

        assert data["hid"] == row_id, (
            f"Expected row ID to match. However, got {data['id']} vs. {row_id}"
        )


@pytest.mark.parametrize("json_option", JSON_OPTIONS)
@pytest.mark.parametrize("list_option", LIST_OPTIONS)
def test_list(config, list_option, json_option):  # noqa: F811
    stdout = _run_cli_command(
        ["data", "list"] + list_option + json_option,
        client=config["client"],
    )

    # check that we can parse into JSON
    if json_option == ["--json"]:
        _check_json(stdout)


@pytest.mark.parametrize("json_option", JSON_OPTIONS)
def test_show_db(config, json_option):  # noqa: F811
    print(config["databases"][0])

    stdout = _run_cli_command(
        ["data", "show", config["databases"][0]] + json_option,
        client=config["client"],
    )

    if json_option == ["--json"]:
        _check_json(stdout)


@beartype
def _check_json(txt: str) -> Union[dict, list[dict]]:
    """helper function to check that a string is valid JSON"""
    data = json.loads(txt)

    return data
