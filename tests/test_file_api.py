"""this module tests the file API"""

import os

import pytest

from deeporigin.platform import file_api
from tests.utils import config  # noqa: F401


@pytest.mark.dependency()
def test_health(config):  # noqa: F811
    """test the health API"""

    data = file_api.check(client=config["client"])
    assert data["status"] == "ok"
    assert data["info"]["file-service"]["status"] == "up"


@pytest.mark.dependency(depends=["test_health"])
def test_get_all_files(config):  # noqa: F811
    """check that there are some files in entities/"""

    files = file_api.get_object_directory(
        file_path="entities/",
        recursive=True,
        client=config["client"],
    )
    assert len(files) > 0, "should be some files in entities/"

    print(f"Found {len(files)} files")


@pytest.mark.dependency(depends=["test_health"])
def test_download_file(config):  # noqa: F811
    """test the file download API"""

    files = file_api.get_object_directory(
        file_path="entities/",
        recursive=True,
        client=config["client"],
    )
    assert len(files) > 0, "should be some files in entities/"

    if config["mock"]:
        return

    # can't run this part of the test in mock mode
    # because file downloads aren't mocked

    local_path = file_api.download_file(
        remote_path=files[0].Key,
        client=config["client"],
    )

    assert os.path.exists(local_path), "should have downloaded the file"
