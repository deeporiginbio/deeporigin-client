"""bridge module to interact with the platform files api"""

import os

from deeporigin.files import FilesClient


def download_files(files: list[str], *, client=None):
    """download files from UFA, if needed, and save to ~/.deeporigin/"""

    if client is None:
        client = FilesClient()

    # Construct src_to_dest mapping for download_files
    home_dir = os.path.expanduser("~")
    src_to_dest = {}
    for remote_path in files:
        local_path = os.path.join(home_dir, ".deeporigin", remote_path)
        os.makedirs(os.path.dirname(local_path), exist_ok=True)
        if not os.path.exists(local_path):
            src_to_dest[remote_path] = local_path

    if src_to_dest:
        client.download_files(src_to_dest)

    # check that every file was downloaded
    for remote_path, local_path in src_to_dest.items():
        if not os.path.exists(local_path):
            raise FileNotFoundError(
                f"File {remote_path} not downloaded to {local_path}"
            )
