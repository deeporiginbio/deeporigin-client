"""bridge module to interact with the platform tools api"""

from concurrent.futures import ThreadPoolExecutor, as_completed
import os
import sys

from beartype import beartype

from deeporigin.platform.utils import _add_functions_to_module

__all__ = _add_functions_to_module(
    module=sys.modules[__name__],
    api_name="FilesApi",
    sdk_name="files",
)


@beartype
def upload_file(*, local_path: str, remote_path: str) -> dict:
    """upload a file using the platform files api"""

    with open(local_path, "rb") as f:
        file_bytes = f.read()

    file_name = os.path.basename(local_path)

    # this function will be imported by _add_functions_to_module
    return put_object(file_path=remote_path, file=(file_name, file_bytes))  # noqa: F821


@beartype
def upload_files(*, files: dict[str, str], max_workers: int = 10):
    """Upload multiple files in parallel using the platform files api.

    Args:
        files: Dict mapping source paths to destination paths
              e.g. {'/path/to/local/file.txt': 'remote/path/file.txt'}
        max_workers: Maximum number of concurrent uploads (default: 10)

    Returns:
        Dict of successfully uploaded files (same format as input)
    """
    successful_uploads = {}

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        # Create future to file mapping
        future_to_file = {
            executor.submit(upload_file, local_path=src, remote_path=dest): (src, dest)
            for src, dest in files.items()
        }

        # Process completed uploads
        for future in as_completed(future_to_file):
            src, dest = future_to_file[future]
            try:
                future.result()
                successful_uploads[src] = dest
            except Exception as e:
                print(f"Error uploading {src}: {str(e)}")

    return successful_uploads
