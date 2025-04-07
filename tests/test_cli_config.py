from deeporigin.config import get_value
from tests.utils import _run_cli_command
from tests.utils import minimal_config as config  # noqa: F401


def test_set_config(config):  # noqa: F811
    org_id = "foo_1234"
    try:
        org_id = get_value()["organization_id"]
    except Exception:
        # this may not work on github actions
        pass

    stdout = _run_cli_command(
        ["config", "set", "organization_id", "foo_3423"],
        client=config["client"],
    )

    assert "foo_3423" in stdout, "Failed to set organization_id"

    # reset to what it was before
    stdout = _run_cli_command(
        ["config", "set", "organization_id", org_id],
        client=config["client"],
    )

    assert org_id in stdout, "Failed to set organization_id"
