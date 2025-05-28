import os

import pytest

from deeporigin.drug_discovery import EXAMPLE_DATA_DIR
from deeporigin.files.files_client import FilesClient
from tests.utils import config  # noqa: F401


def test_upload_and_download_files(config, tmp_path):  # noqa: F811
    if config["mock"]:
        pytest.skip("test skipped with mock client")

    files_client = FilesClient()

    files = os.listdir(EXAMPLE_DATA_DIR)
    src_to_dest = {
        os.path.join(EXAMPLE_DATA_DIR, file): os.path.join("test-upload/", file)
        for file in files
    }

    success, _ = files_client.upload_files(src_to_dest)
    assert success, "Failed to upload files"

    # Prepare download mapping: remote dest -> local temp dir
    download_dest = {
        remote: str(tmp_path / os.path.basename(remote))
        for remote in src_to_dest.values()
    }
    success, _ = files_client.download_files(download_dest)
    assert success, "Failed to download files"

    for file in files:
        assert os.path.exists(tmp_path / file), (
            f"File {file} not found in temp download directory"
        )


def test_list_folder(config):  # noqa: F811
    if config["mock"]:
        pytest.skip("test skipped with mock client")

    files_client = FilesClient()

    files = os.listdir(EXAMPLE_DATA_DIR)
    src_to_dest = {
        os.path.join(EXAMPLE_DATA_DIR, file): os.path.join("test-upload/", file)
        for file in files
    }

    success, _ = files_client.upload_files(src_to_dest)
    assert success, "Failed to upload files"

    data = files_client.list_folder("test-upload/")

    assert set(data.keys()) == set(src_to_dest.values()), (
        "Failed to list files correctly"
    )


def test_upload_files_missing_local_files(config):  # noqa: F811
    if config["mock"]:
        pytest.skip("test skipped with mock client")

    files_client = FilesClient()

    # Use a nonexistent file path
    missing_file = "/tmp/this_file_does_not_exist_123456789.txt"
    src_to_dest = {missing_file: "test-upload/missing.txt"}
    with pytest.raises(FileNotFoundError) as excinfo:
        files_client.upload_files(src_to_dest)

    assert missing_file in str(excinfo.value)
