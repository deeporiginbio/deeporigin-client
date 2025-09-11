"""this module tests the file API"""

import os

import pytest

from deeporigin.platform import Client, file_api
from tests.utils import config  # noqa: F401


@pytest.mark.dependency()
def test_health(config):  # noqa: F811
    """test the health API"""

    if config["mock"]:
        pytest.skip("test skipped with mock client")

    data = file_api.check(client=Client())
    assert data["status"] == "ok"
    assert data["info"]["file-service"]["status"] == "up"


@pytest.mark.dependency(depends=["test_health"])
def test_get_all_files(config):  # noqa: F811
    """check that there are some files in entities/"""

    if config["mock"]:
        pytest.skip("test skipped with mock client")

    files = file_api.get_object_directory(file_path="entities/", recursive=True)
    assert len(files) > 0, "should be some files in entities/"

    print(f"Found {len(files)} files")


@pytest.mark.dependency(depends=["test_health"])
def test_download_file(config):  # noqa: F811
    """test the file download API"""

    if config["mock"]:
        pytest.skip("test skipped with mock client")

    files = file_api.get_object_directory(file_path="entities/", recursive=True)
    assert len(files) > 0, "should be some files in entities/"

    local_path = file_api.download_file(remote_path=files[0].Key)

    assert os.path.exists(local_path), "should have downloaded the file"
