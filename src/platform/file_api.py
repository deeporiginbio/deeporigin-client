"""bridge module to interact with the platform files api"""

import concurrent.futures
import os
import sys
from typing import Optional

from beartype import beartype

from deeporigin.platform.utils import _add_functions_to_module
from deeporigin.utils.core import _ensure_do_folder

__all__ = _add_functions_to_module(
    module=sys.modules[__name__],
    api_name="FileApi",
)


DO_FOLDER = _ensure_do_folder()


@beartype
def upload_file(*, local_path: str, remote_path: str, client=None):
    """upload a single file to UFA

    Args:
        local_path (str): The local path of the file to upload.
        remote_path (str): The remote path of the file to upload.
    """

    with open(local_path, "rb") as f:
        file_data = (remote_path, f.read())
        put_object(file_path=remote_path, file=file_data, client=client)  # noqa: F821


@beartype
def download_file(
    *,
    remote_path: str,
    local_path: Optional[str] = None,
    client=None,
):
    """download a single file from UFA to ~/.deeporigin/, or some other local path

    Args:
        remote_path (str): The remote path of the file to download.
        local_path (str): The local path to save the file to. If None, uses ~/.deeporigin/.
    """

    if local_path is None:
        local_path = os.path.join(DO_FOLDER, remote_path)

    os.makedirs(os.path.dirname(local_path), exist_ok=True)

    response = get_signed_url(  # noqa: F821
        file_path=remote_path,
        client=client,
    )

    from deeporigin.utils.network import download_sync

    download_sync(response.url, save_path=local_path)


@beartype
def upload_files(
    files: dict[str, str],
    *,
    client=None,
):
    """Upload multiple files in parallel. files: {local_path: remote_path}. Raises RuntimeError if any upload fails.

    Args:
        files (dict[str, str]): A dictionary of local paths to remote paths.
    """

    results = []
    errors = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        future_to_pair = {
            executor.submit(
                upload_file, local_path=lp, remote_path=rp, client=client
            ): (lp, rp)
            for lp, rp in files.items()
        }
        for future in concurrent.futures.as_completed(future_to_pair):
            lp, rp = future_to_pair[future]
            try:
                results.append(future.result())
            except Exception as e:
                errors.append((lp, rp, e))
    if errors:
        error_msgs = "\n".join(
            [
                f"Upload failed for local_path={lp}, remote_path={rp}: {str(err)}"
                for lp, rp, err in errors
            ]
        )
        raise RuntimeError(f"Some uploads failed in upload_files:\n{error_msgs}")
    return results


@beartype
def download_files(files: dict[str, str | None], *, client=None):
    """Download multiple files in parallel. files: {remote_path: local_path or None}. If local_path is None, use default. Raises RuntimeError if any download fails."""
    results = []
    errors = []
    with concurrent.futures.ThreadPoolExecutor() as executor:
        future_to_pair = {
            executor.submit(
                download_file, remote_path=rp, local_path=lp, client=client
            ): (rp, lp)
            for rp, lp in files.items()
        }
        for future in concurrent.futures.as_completed(future_to_pair):
            rp, lp = future_to_pair[future]
            try:
                results.append(future.result())
            except Exception as e:
                errors.append((rp, lp, e))
    if errors:
        error_msgs = "\n".join(
            [
                f"Download failed for remote_path={rp}, local_path={lp}: {str(err)}"
                for rp, lp, err in errors
            ]
        )
        raise RuntimeError(f"Some downloads failed in download_files:\n{error_msgs}")
    return results
