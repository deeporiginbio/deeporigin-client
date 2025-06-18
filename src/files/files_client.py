"""
DeepOrigin FilesClient Module

This module provides a client for interacting with the DeepOrigin file service.
It offers high-level methods for file operations such as upload, download,
listing and synchronizing directories between local and remote storage.

Example usage:
    ```python
    from deeporigin.files import FilesClient

    # Create a client (uses default authentication if token not provided)
    client = FilesClient()

    # List files in a directory
    files = client.list_folder('data/')
    print(f"Found {len(files)} files")

    # Download a file from remote storage to local path
    success = client.download_file(
        src="path/to/remote/file.txt",
        dest="local/path/file.txt"
    )

    if success:
        print("File downloaded successfully")
    else:
        print("Download failed")
    ```
"""

import concurrent.futures
from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum

# import io
import logging
import os
import time
from typing import Callable

from deeporigin import auth
from deeporigin.utils.config import get_value
from deeporigin.utils.network import _get_domain_name

from .file_service import generated_api as genapi
from .file_service.schema_api_httpx_client import FilePayload, SchemaApiHttpxClient

# Remote path prefix for file service URLs
# REMOTE_PATH_PREFIX = "files:///"
# FILES_URL_PATTERN = re.compile(rf"^{REMOTE_PATH_PREFIX}(?P<path>.*)$")
REMOTE_PATH_PREFIX = ""
# The intent of path-prefix was to distinguish between local and remote paths.
# This is only really useful if multiple destinations were supported by the same
# API, if this happens we may re-enable this. The down side is that it complicates
# use; how if we introduce it may (a) report it in all paths as part of list results
# and/or (b) consider setting it or the whole Client instead of every url.


class FolderSyncMode(Enum):
    """
    Enum for folder synchronization modes.
    """

    # We don't know if the destination file not present locally was deleted remotely or not.
    # This means we can't do "proper" sync, we just have different overwrite modes.

    REPLACE = "replace_older"  # replace older files, delete files not present in source
    MERGE = "merge"  # replace older files, does NOT delete files not present in source
    FORCE_OVERWRITE = (
        "force_overwrite"  # replace ALL files, delete files not present in source
    )
    COPY_NON_EXISTING = (
        "copy_non_existing"  # copy files that don't exist in destination
    )

    @classmethod
    def default(cls) -> "FolderSyncMode":
        """Return the default sync mode (REPLACE)"""
        return cls.REPLACE

    def can_delete_nonexisting(self) -> bool:
        return self == FolderSyncMode.REPLACE or self == FolderSyncMode.FORCE_OVERWRITE


@dataclass
class FileMetadata:
    """
    Structured representation of file metadata, similar to S3 metadata.

    This class provides a typed interface to file metadata while maintaining
    compatibility with dictionary-based access.
    """

    # Common S3-like metadata fields
    key_path: str | None = None  # File path or its key in the storage system
    last_modified: str | None = None
    etag: str | None = None
    size: int | None = None  # Size of file/object in bytes
    storage_class: str | None = None
    content_type: str | None = None

    # Store the original dictionary for access to all fields
    raw_dict: dict[str, any] = field(default_factory=dict)

    # Mapping between S3/HTTP header names and class attribute names - shared across instances
    # Note: Metadata and list dictionary return different fields, so we have two mappings.
    # For our request content-length is size, but technically may be size of HTTP instead
    # if is content, so may need to adjust with protocol in the future.
    METADATA_HEAD_FIELD_MAPPING = {
        "key": "key_path",  # Only comes in list (not HEAD)
        "last-modified": "last_modified",
        "x-last-modified": "last_modified",  # comes in different in list and HEAD
        "etag": "etag",
        "content-length": "size",
        "content-type": "content_type",
        "x-amz-storage-class": "storage_class",
    }
    METADATA_LIST_FIELD_MAPPING = {
        "key": "key_path",  # Only comes in list (not HEAD)
        "etag": "etag",
        "lastmodified": "last_modified",  # LastModified -> last_modified
        "size": "size",  # Size -> size  (in list)
        "contenttype": "content_type",  # ContentType -> content_type
        "storageclass": "storage_class",  # StorageClass -> storage_c
    }

    @property
    def is_folder(self) -> bool:
        """
        Determine if this object represents a folder.

        Returns:
            True if the object appears to be a folder based on:
            - Path ending with '/'
            - Content-type being 'application/x-directory'
        """
        if self.key_path and self.key_path.endswith("/"):
            return True
        if self.content_type == "application/x-directory":
            return True
        return False

    @classmethod
    def from_dict(
        cls,
        data: dict[str, any],
        *,
        key_path: str | None = None,
        mapping: dict[str, str] = METADATA_HEAD_FIELD_MAPPING,
    ) -> "FileMetadata":
        """
        Create a FileMetadata instance from a dictionary.
        Args:
            data: Dictionary containing metadata fields
        Returns:
            FileMetadata instance with fields populated from the dictionary
        """
        instance = cls()
        instance.raw_dict = data  # Store reference, no need to copy

        # Process each key in the input dictionary, setting matching attributes
        for key, value in data.items():
            key_lower = key.lower()
            attr_name = mapping.get(key_lower, key)

            if hasattr(instance, attr_name):
                if attr_name == "size":
                    try:
                        value = int(value)
                    except (ValueError, TypeError):
                        pass
                setattr(instance, attr_name, value)

        if key_path:
            instance.key_path = key_path
        return instance

    def to_updated_dict(
        self, *, mapping: dict[str, str] = METADATA_HEAD_FIELD_MAPPING
    ) -> dict[str, any]:
        """
        Get a copy of the original dictionary with updated values from attributes.
        Maintains original key capitalization.
        Returns:
            Updated copy of the original metadata dictionary
        """
        result = self.raw_dict.copy()

        # This assumes _FIELD_MAPPING includes all needed attributes.
        # For each _FIELD_MAPPING entry, if the attribute exists assign it as string.
        # TBD -- may result in repeated keys if there is repeated field mapping
        #   (ex: last-midified)
        for header_key, attr_name in self._FIELD_MAPPING.items():
            if hasattr(self, attr_name) and getattr(self, attr_name) is not None:
                result[header_key] = str(getattr(self, attr_name))
        return result

    def get_last_modified_timestamp(self) -> float | None:
        """
        Get the last_modified field as a timestamp.
        Returns:
            Timestamp as float, or None if last_modified is not available or invalid
        """
        if not self.last_modified:
            return None
        try:
            # Parse ISO format timestamp
            dt = datetime.fromisoformat(self.last_modified.replace("Z", "+00:00"))
            return dt.timestamp()
        except (ValueError, TypeError):
            return None


class FilesClient:
    """
    Client for interacting with the DeepOrigin file service.

    This client provides high-level methods for file operations like upload, download,
    listing folder, and synchronizing files between local and remote storage.
    """

    # FilesClient instance attributes
    base_url: str
    organization_id: str
    client: SchemaApiHttpxClient
    # Logger used for client logging. Re-assign if you need to override.
    # We'll generally log errors if expected i/o fails, warnings on retries
    logger: logging.Logger
    # Function to be called to create a progress callback for multi-file i/o operations.
    # None by default, tqdm for a Jupyter notebook. To turn off set to None after constructor.
    progress_callback_create_func: (
        Callable[[str], Callable[[int, int], None] | None] | None
    )
    # These control retries for upload and download
    io_max_retries: int
    io_retry_delay: float

    def __init__(
        self,
        base_url: str | None = None,
        token: str | None = None,
        verify_ssl: bool | str = True,
        organization_id: str | None = None,
    ):
        """
        Initialize the FilesClient.
        Args:
            base_url: Base URL of the file service API
            token: Authentication token (optional)
            verify_ssl: Whether to verify SSL certificates or path to CA bundle
        """

        self.logger = logging.getLogger(__name__)

        self.io_max_retries = 3
        self.io_retry_delay = 0.5  # retry after this many seconds

        self.progress_callback_create_func = None
        if _is_notebook():
            self.progress_callback_create_func = create_files_tqdm_callback

        if not organization_id:
            self.organization_id = get_value()["organization_id"]
        else:
            self.organization_id = organization_id

        if not base_url:
            base_url = _get_domain_name()
        self.base_url = base_url
        if not token:
            tokens = auth.get_tokens(refresh=False)
            token = tokens["access"]

        # TBD: Do we want to still check null token and maybe throw exception?
        #      self.client = AuthenticatedClient(
        #     base_url=self.base_url, token=token, verify_ssl=verify_ssl
        # )

        schema_path = os.path.join(
            os.path.dirname(os.path.abspath(__file__)),
            "file_service/gen_file_openapi.yaml",
        )
        # Initialize the simple API client
        self.client = SchemaApiHttpxClient(
            base_url=self.base_url,
            token=token,
            verify_ssl=verify_ssl,
            schema_path=schema_path,
        )

    def _extract_path_from_url(self, url: str) -> str:
        # This was used in the past to extract the path from the URL, seems unnecessary.
        # match = FILES_URL_PATTERN.match(url)
        # return match.group("path")
        return url

    def list_folder(
        self, path: str, recursive: bool = False
    ) -> dict[str, FileMetadata]:
        """
        List files and folders at the specified path. Returned dictionary keys will have
        the full remote path, so the resulting path key can be used directly for download, etc.
        Args:
            path: Remote path where to list files from
            recursive: Optional "recursive" flag. If specified will include sub-folder content.
        Returns:
            Dictionary mapping file paths to FileMetadata objects, maintaining the order
            from the underlying list response. Key string will include the argument
            path part; it is not removed.
        """
        remote_path = self._extract_path_from_url(path)
        result = {}
        continuation_token = None
        max_iterations = 10
        current_iteration = 0

        while current_iteration < max_iterations:
            call_args = genapi._get_get_object_kwargs(
                orgFriendlyId=self.organization_id,
                filePath=remote_path,
                list_type=2,
                recursive="true" if recursive else None,
                continuation_token=continuation_token,
            )
            response = self.client.httpx_request_and_process(call_args)

            if response.status_code == 200:
                raw_files = self._parse_list_response(response.content)
                raw_files_data = raw_files["data"]

                # Add files from this batch to the result
                for file_data in raw_files_data:
                    metadata = FileMetadata.from_dict(
                        file_data, mapping=FileMetadata.METADATA_LIST_FIELD_MAPPING
                    )
                    if metadata.key_path:
                        result[metadata.key_path] = metadata

                # Check for continuation token
                continuation_token = raw_files.get("continuationToken")
                if not continuation_token:
                    break  # No more pages to fetch

                current_iteration += 1
            else:
                self.logger.error(
                    f"Failed to list directory {path}: {response.status_code}"
                )
                break

        if current_iteration >= max_iterations:
            self.logger.warning(
                f"Reached maximum iteration limit ({max_iterations}) while listing directory {path}"
            )

        return result

    def _parse_list_response(self, content: bytes) -> list[dict[str, any]]:
        """
        Parse the response from a list directory operation.
        Args:
            content: Response content bytes
        Returns:
            List of file/directory metadata objects
        """
        # This implementation will depend on the actual response format
        # For now, we'll return a placeholder
        import json

        try:
            return json.loads(content)
        except json.JSONDecodeError:
            self.logger.error("Could not parse list response as JSON")
            return []

    def _remote_file_exists(self, path: str) -> bool:
        # If we got metadata without an error, the file exists
        try:
            metadata = self.get_metadata(path)
            return bool(metadata)
        except Exception:
            # If there was an error getting metadata, the file probably doesn't exist
            return False

    def create_folder(self, path: str) -> bool:
        """
        Create a folder in the remote storage.
        In case of S3, should create a zero-byte object with a trailing slash and special
        content type to represent a folder in object storage.
        Args:
            path: Remote path where to create the folder. Must end with '/'
        Returns:
            True if folder creation was successful, False otherwise
        """
        if not path.endswith("/"):
            error_msg = f"Folder path must end with '/': {path}"
            self.logger.error(error_msg)
            raise ValueError(error_msg)

        remote_path = self._extract_path_from_url(path)

        try:
            # TBD: ideally mime-type created would be "application/x-directory"
            call_args = genapi._get_put_object_kwargs(
                orgFriendlyId=self.organization_id,
                filePath=remote_path,
                # file=file_obj,
            )
            response = self.client.httpx_request_and_process(call_args)

            if response.status_code == 200:
                return True
            else:
                self.logger.error(
                    f"Failed to create folder {path}: Status code {response.status_code}"
                )
                return False
        except Exception as e:
            self.logger.error(f"Failed to create folder {path}: {e}")
            return False

    def delete_folder(self, path: str) -> bool:
        """
        Delete a folder from remote storage. This is not recursive and will only
        delete the folder itself, not its contents.
        Args:
            path: Path to the folder. Must end with '/'
        Returns:
            True if deletion was successful, False otherwise
            ValueError: If path does not end with '/'
        """
        if not path.endswith("/"):
            error_msg = f"Folder path must end with '/': {path}"
            self.logger.error(error_msg)
            raise ValueError(error_msg)

        # Use the existing delete_file method as the operation is the same
        return self._delete_single_file(path)[0]

    def upload_file(self, src: str, dest: str) -> bool:
        """
        Upload a file from local path to remote storage. File is overwritten if it exists.
        Args:
            src: Local source path
            dest: Remote destination path
        Returns:
            True if upload was successful, False otherwise
        """
        # Use the new upload_files method with a single file
        success, _ = self.upload_files({src: dest})
        return success

    def upload_files(
        self,
        src_to_dest: dict[str, str],
        progress_callback: Callable[[int, int], None] | None = None,
    ) -> tuple[bool, dict[str, str]]:
        """
        Upload multiple files from local paths to remote storage paths. Upload is
        done in parallel and may execute out of order.
        Args:
            src_to_dest: Dictionary mapping source file paths to destination paths
        Returns:
            Tuple containing:
                - Boolean indicating if all uploads were successful
                - Dictionary mapping source paths to status ("OK" or error message)
        """

        missing_files = [src for src in src_to_dest.keys() if not os.path.exists(src)]
        if missing_files:
            raise FileNotFoundError(
                f"The following files were not found: {', '.join(missing_files)}"
            )

        if progress_callback is None and self.progress_callback_create_func:
            progress_callback = self.progress_callback_create_func("Uploading")

        return self._process_files_in_parallel(
            src_to_dest,
            self._upload_single_file,
            "Upload",
            progress_callback,
        )

    def download_files(
        self,
        src_to_dest: dict[str, str],
        progress_callback: Callable[[int, int], None] | None = None,
    ) -> tuple[bool, dict[str, str]]:
        """
        Download multiple files from remote storage to local paths. Download is
        done in parallel and may execute out of order.
        Args:
            src_to_dest: Dictionary mapping source file paths to destination paths
        Returns:
            Tuple containing:
                - Boolean indicating if all downloads were successful
                - Dictionary mapping source paths to status ("OK" or error message)
        """
        if progress_callback is None and self.progress_callback_create_func:
            progress_callback = self.progress_callback_create_func("Downloading")

        return self._process_files_in_parallel(
            src_to_dest, self._download_single_file, "Download", progress_callback
        )

    def download_file(self, src: str, dest: str) -> bool:
        """
        Download a file from remote storage to local path.
        Args:
            src: Remote source path
            dest: Local destination path
        Returns:
            True if download was successful, False otherwise
        """
        success, _ = self.download_files({src: dest})
        return success

    def delete_file(self, path: str) -> bool:
        """
        Delete a file from remote storage.
        Args:
            path: Path in the format
        Returns:
            True if deletion was successful, False otherwise
        """
        success, _ = self.delete_files([path])
        return success

    def delete_files(
        self,
        paths: list[str],
        progress_callback: Callable[[int, int], None] | None = None,
    ) -> tuple[bool, dict[str, str]]:
        """
        Delete multiple files from remote storage in parallel.
        Args:
            paths: List of file paths to delete
        Returns:
            Tuple containing:
                - Boolean indicating if all deletions were successful
                - Dictionary mapping paths to status ("DELETED" or error message)
        """
        if progress_callback is None and self.progress_callback_create_func:
            progress_callback = self.progress_callback_create_func("Deleting")

        return self._process_files_in_parallel(
            paths, self._delete_single_file, "Delete", progress_callback
        )

    def _process_files_in_parallel(
        self,
        _args_list: dict[str, str] | list[str],
        single_file_io_func,
        operation_name: str,
        progress_callback: Callable[[int, int], None] | None = None,
    ) -> tuple[bool, dict[str, str]]:
        if not _args_list:
            return True, {}

        results = {}
        all_success = True
        total_items = len(_args_list)  # Get total count
        items_completed = 0  # Initialize counter

        # Determine the iterable based on input type
        args_list = _args_list.items() if isinstance(_args_list, dict) else _args_list

        # For a single file, use the function directly (no loop, callback called once)
        if total_items == 1:
            src_dest_args = next(iter(args_list))
            src = (
                src_dest_args[0] if isinstance(src_dest_args, tuple) else src_dest_args
            )
            success, status = single_file_io_func(src_dest_args)
            items_completed = 1
            if progress_callback:
                progress_callback(items_completed, total_items)
            results[src] = status
            return success, results
        else:
            # For multiple files, use concurrent futures
            with concurrent.futures.ThreadPoolExecutor() as executor:
                future_to_src = {}

                for src_dest_args in args_list:
                    future = executor.submit(single_file_io_func, src_dest_args)
                    src = (
                        src_dest_args[0]
                        if isinstance(src_dest_args, tuple)
                        else src_dest_args
                    )
                    future_to_src[future] = src

                for future in concurrent.futures.as_completed(future_to_src):
                    src = future_to_src[future]
                    try:
                        success, status = future.result()
                        results[src] = status
                        if not success:
                            all_success = False

                    except Exception as e:
                        results[src] = f"{operation_name} failed with error: {str(e)}"
                        self.logger.error(results[src])
                        all_success = False
                    finally:
                        # Increment counter and call the callback
                        items_completed += 1
                        if progress_callback:
                            progress_callback(items_completed, total_items)

        return all_success, results

    def _upload_single_file(
        self, str_to_dest_tuple: tuple[str, str]
    ) -> tuple[bool, str]:
        src, dest = str_to_dest_tuple
        if not os.path.exists(src):
            return False, f"Source file not found: {src}"

        remote_path = self._extract_path_from_url(dest)
        last_error = None

        for attempt in range(self.io_max_retries):
            try:
                with open(src, "rb") as f:
                    file_obj = FilePayload(
                        payload=f,
                        file_name=os.path.basename(src),
                        mime_type="application/octet-stream",
                    )
                    # format of expected dictonary from file_obj is
                    # {"file": (file_name, file_content, mime_type) }

                    call_args = genapi._get_put_object_kwargs(
                        orgFriendlyId=self.organization_id,
                        filePath=remote_path,
                        file=file_obj,
                    )
                    response = self.client.httpx_request_and_process(call_args)

                if response.status_code == 200:
                    return True, "OK"

                # Determine if we should retry based on status code. 500 codes are server errors.
                # 429 is rate limit exceeded, 408 is request timeout.
                if response.status_code in [429, 408] or (
                    500 <= response.status_code < 600
                ):
                    if attempt < self.io_max_retries - 1:
                        last_error = f"Attempt {attempt + 1} failed with status {response.status_code}"
                        self.logger.warning(
                            f"{last_error}, retrying {src} upload in {self.io_retry_delay}s..."
                        )
                        time.sleep(self.io_retry_delay)
                        continue

                # For other status codes, don't retry
                last_error = f"Upload failed with status code: {response.status_code}"
                break

            except (ConnectionError, TimeoutError) as e:
                # Retry on connection errors
                if attempt < self.io_max_retries - 1:
                    last_error = f"Connection error on attempt {attempt + 1}: {str(e)}"
                    self.logger.warning(
                        f"{last_error}, retrying {src} upload in {self.io_retry_delay}s"
                    )
                    time.sleep(self.io_retry_delay)
                    continue
                last_error = f"Upload failed with connection error: {str(e)}"
                break

            except Exception as e:
                last_error = f"Upload failed with error: {str(e)}"
                break

        # We've exhausted all retries, return the last error message
        self.logger.error(last_error)
        return False, last_error

    def _download_single_file(
        self, str_to_dest_tuple: tuple[str, str]
    ) -> tuple[bool, str]:
        src, dest = str_to_dest_tuple
        remote_path = self._extract_path_from_url(src)
        last_error = None

        for attempt in range(self.io_max_retries):
            try:
                # TBD: The get_object API will need to be called differently for file download
                # vs. directory listing. This implementation assumes that omitting list_type
                # will result in a file download.

                call_args = genapi._get_get_object_kwargs(
                    orgFriendlyId=self.organization_id, filePath=remote_path
                )
                response = self.client.httpx_request_and_process(call_args)
                # list_type=0.0  - wasn't working # Assuming 0.0 means download file
                # looks like bad if(type) check in file-service controller --TBD Niels

                if response.status_code == 200:
                    # Ensure the directory exists
                    os.makedirs(os.path.dirname(os.path.abspath(dest)), exist_ok=True)

                    # Write the content to the destination file
                    with open(dest, "wb") as f:
                        f.write(response.content)
                    return True, "OK"

                    # Determine if we should retry based on status code. 500 codes are server errors.
                    # 429 is rate limit exceeded, 408 is request timeout.
                    if response.status_code in [429, 408] or (
                        500 <= response.status_code < 600
                    ):
                        if attempt < self.io_max_retries - 1:
                            last_error = f"Attempt {attempt + 1} failed with status {response.status_code}"
                            self.logger.warning(
                                f"{last_error}, retrying {src} download in {self.io_retry_delay}s..."
                            )
                            time.sleep(self.io_retry_delay)
                            continue

                    # For other status codes, don't retry
                    last_error = (
                        f"Download failed with status code: {response.status_code}"
                    )
                    break

            except (ConnectionError, TimeoutError) as e:
                # Retry on connection errors
                if attempt < self.io_max_retries - 1:
                    last_error = f"Connection error on attempt {attempt + 1}: {str(e)}"
                    self.logger.warning(
                        f"{last_error}, retrying {src} download in {self.io_retry_delay}s"
                    )
                    time.sleep(self.io_retry_delay)
                    continue
                last_error = f"Download failed with connection error: {str(e)}"
                break

            except Exception as e:
                last_error = f"Download failed with error: {str(e)}"
                break

        # We've exhausted all retries, return the last error message
        self.logger.error(last_error)
        return False, last_error

    def _delete_single_file(self, path: str) -> tuple[bool, str]:
        remote_path = self._extract_path_from_url(path)
        last_error = None

        for attempt in range(self.io_max_retries):
            try:
                call_args = genapi._get_delete_object_kwargs(
                    orgFriendlyId=self.organization_id, filePath=remote_path
                )
                response = self.client.httpx_request_and_process(call_args)

                # Successful deletion (200 OK or 204 No Content)
                # Also treat 404 Not Found as success (idempotent delete)
                if response.status_code in [200, 204]:
                    return True, "DELETED"

                # Determine if we should retry based on status code
                if response.status_code in [429, 408] or (
                    500 <= response.status_code < 600
                ):
                    if attempt < self.io_max_retries - 1:
                        last_error = f"Attempt {attempt + 1} failed with status {response.status_code}"
                        self.logger.warning(
                            f"{last_error}, retrying {path} deletion in {self.io_retry_delay}s..."
                        )
                        time.sleep(self.io_retry_delay)
                        continue
                # For other status codes, don't retry
                last_error = f"Delete failed with status code: {response.status_code}"
                break

            except (ConnectionError, TimeoutError) as e:
                # Retry on connection errors
                if attempt < self.io_max_retries - 1:
                    last_error = f"Connection error on attempt {attempt + 1}: {str(e)}"
                    self.logger.warning(
                        f"{last_error}, retrying {path} deletion in {self.io_retry_delay}s"
                    )
                    time.sleep(self.io_retry_delay)
                    continue
                last_error = f"Delete failed with connection error: {str(e)}"
                break

            except Exception as e:
                last_error = f"Delete failed with error: {str(e)}"
                break

        # We've exhausted all retries, return the last error message
        self.logger.error(last_error)
        return False, last_error

    def get_metadata(self, path: str) -> FileMetadata | None:
        """
        Get metadata about a file or directory.
        Args:
            path: Path on the remote server
        Returns:
            FileMetadata object containing the metadata
        """
        remote_path = self._extract_path_from_url(path)

        call_args = genapi._get_head_object_kwargs(
            orgFriendlyId=self.organization_id, filePath=remote_path
        )
        response = self.client.httpx_request_and_process(call_args)

        if response.status_code == 200:
            # Extract metadata from response headers and create a FileMetadata object
            return FileMetadata.from_dict(dict(response.headers), key_path=remote_path)
        else:
            self.logger.error(
                f"Failed to get metadata for {path}: {response.status_code}"
            )
            return None

    def sync_folder_down(
        self,
        remote_path: str,
        local_path: str,
        mode: FolderSyncMode = FolderSyncMode.REPLACE,
        progress_callback: Callable[[int, int], None] | None = None,
    ) -> tuple[bool, dict[str, str]]:
        """
        Sync files from remote to local directory.
        Args:
            remote_path: Remote path on the server
            local_path: Local directory path
            mode: Synchronization mode; REPLACE by default, indicating that older files are
            overwritten by newer remote files, and non-existing remote files are deleted locally.
        Returns:
            Tuple containing:
                - Boolean True if the sync of all files was successful, False if there was any error.
                  Need to look at the dictonary entry to see what happened per file.
                - Dictionary mapping file paths to their status. Possibles success status include:
                  "OK", "DELETED", "KEPT". Othewise the string will contain the error message.
        """
        # Get list of remote files with metadata
        # todo: add/figure out recursive flag
        remote_files_dict = self.list_folder(remote_path, recursive=True)

        # Extract the base remote path without the prefix for path comparisons
        base_remote_path = self._extract_path_from_url(remote_path)

        # Ensure local directory exists
        os.makedirs(local_path, exist_ok=True)

        # Track all file operations
        file_statuses = {}
        success = True

        # Build a dictionary of files to download
        files_to_download = {}

        # Process remote files
        for remote_file_key, file_metadata in remote_files_dict.items():
            # Construct paths
            full_remote_file_path = f"{REMOTE_PATH_PREFIX}{remote_file_key}"
            remote_file_key_tail = os.path.relpath(remote_file_key, base_remote_path)
            local_file_path = os.path.join(local_path, remote_file_key_tail)

            # Ensure local directory exists
            local_dir = os.path.dirname(local_file_path)
            if local_dir:
                os.makedirs(local_dir, exist_ok=True)

            # Check if local file exists
            if os.path.exists(local_file_path):
                local_mtime = os.path.getmtime(local_file_path)
                remote_mtime = file_metadata.get_last_modified_timestamp() or 0

                # Determine if we should update the file based on mode
                should_update = False
                if mode == FolderSyncMode.FORCE_OVERWRITE:
                    should_update = True
                elif mode == FolderSyncMode.REPLACE:
                    should_update = local_mtime < remote_mtime
                elif mode == FolderSyncMode.MERGE:
                    should_update = local_mtime < remote_mtime
                elif mode == FolderSyncMode.COPY_NON_EXISTING:
                    should_update = False

                if should_update:
                    # Add to files to download
                    files_to_download[full_remote_file_path] = local_file_path
                else:
                    file_statuses[local_file_path] = "KEPT"
            else:
                # File doesn't exist locally, add to files to download
                files_to_download[full_remote_file_path] = local_file_path

        # Download files in parallel
        if files_to_download:
            download_success, download_results = self.download_files(
                files_to_download, progress_callback=progress_callback
            )
            if not download_success:
                success = False

            # Add download results to file statuses
            for src, status in download_results.items():
                dest = files_to_download[src]
                file_statuses[dest] = status

        # Handle deletion of files in destination that don't exist in source
        if mode.can_delete_nonexisting():
            # Get all local files
            local_files = []
            for root, _, files in os.walk(local_path):
                for filename in files:
                    local_file_path = os.path.join(root, filename)
                    rel_path = os.path.relpath(local_file_path, local_path)
                    local_files.append(rel_path)

            # Find files to delete
            for local_file in local_files:
                local_file_in_remote = os.path.join(base_remote_path, local_file)
                if local_file_in_remote not in remote_files_dict:
                    full_local_path = os.path.join(local_path, local_file)
                    try:
                        os.remove(full_local_path)
                        file_statuses[full_local_path] = "DELETED"
                    except Exception as e:
                        self.logger.error(
                            f"Failed to delete file {full_local_path}: {e}"
                        )
                        file_statuses[full_local_path] = (
                            f"Error deleting file: {str(e)}"
                        )
                        success = False

        return success, file_statuses

    def sync_folder_up(
        self,
        local_path: str,
        remote_path: str,
        mode: FolderSyncMode = FolderSyncMode.REPLACE,
        progress_callback: Callable[[int, int], None] | None = None,
    ) -> tuple[bool, dict[str, str]]:
        """
        Sync files from local to remote directory.
        Args:
            local_path: Local directory path
            remote_path: Remote path on server
            mode: Synchronization mode; REPLACE by default, indicating that older remote files are
            overwritten by newer local files, and non-existing local files are deleted on remote.
        Returns:
            Tuple containing:
                - Boolean True if the sync of all files was successful, False if there was any error.
                  Need to look at the dictonary entry to see what happened per file.
                - Dictionary mapping file paths to their status. Possibles success status include:
                  "OK", "DELETED", "KEPT". Otherwise the string will contain the error message.
        """
        # Extract the base remote path without the prefix for path comparisons
        base_remote_path = self._extract_path_from_url(remote_path)

        # Get list of remote files with metadata to compare
        remote_files_dict = {}
        try:
            remote_files_dict = self.list_folder(
                remote_path
            )  # todo: add recursive flag
        except Exception as e:
            # MA: TBD - the above won't actually re-throw exception
            # if we care bount this -- need special list function or parameter
            self.logger.warning(
                f"Could not list remote directory, assuming it's empty: {e}"
            )

        # Track all file operations
        file_statuses = {}
        success = True

        # Build a dictionary of files to upload
        files_to_upload = {}
        local_rel_paths = set()  # Initialize set to store relative paths

        # Process local files
        for root, _, files in os.walk(local_path):
            for filename in files:
                local_file_path = os.path.join(root, filename)

                # Get relative path from the local_path
                rel_path = os.path.relpath(local_file_path, local_path)
                local_rel_paths.add(rel_path)  # Add for deletion loop below

                remote_file_path = f"{remote_path}/{rel_path}"
                rel_key_path = base_remote_path + "/" + rel_path

                # Check if remote file exists
                if rel_key_path in remote_files_dict:
                    local_mtime = os.path.getmtime(local_file_path)
                    remote_mtime = (
                        remote_files_dict[rel_key_path].get_last_modified_timestamp()
                        or 0
                    )

                    # Determine if we should update the file based on mode
                    should_update = False
                    if mode == FolderSyncMode.FORCE_OVERWRITE:
                        should_update = True
                    elif mode == FolderSyncMode.REPLACE:
                        should_update = local_mtime > remote_mtime
                    elif mode == FolderSyncMode.MERGE:
                        should_update = local_mtime > remote_mtime
                    elif mode == FolderSyncMode.COPY_NON_EXISTING:
                        should_update = False

                    if should_update:
                        # Add to files to upload
                        files_to_upload[local_file_path] = remote_file_path
                    else:
                        file_statuses[remote_file_path] = "KEPT"
                else:
                    # File doesn't exist remotely, add to files to upload
                    files_to_upload[local_file_path] = remote_file_path

        # Upload files in parallel
        if files_to_upload:
            upload_success, upload_results = self.upload_files(
                files_to_upload, progress_callback=progress_callback
            )
            if not upload_success:
                success = False

            # Add upload results to file statuses
            for src, status in upload_results.items():
                dest = files_to_upload[src]
                file_statuses[dest] = status

        # Handle deletion of files in destination that don't exist in source
        if mode.can_delete_nonexisting():
            # local_rel_paths is already populated from the first loop
            # Find remote files to delete
            files_to_delete_list = []
            for remote_file_key in remote_files_dict:
                remote_rel_path = os.path.relpath(remote_file_key, base_remote_path)
                if remote_rel_path not in local_rel_paths:
                    # Construct the full remote path using the original input format
                    full_remote_path = f"{remote_path}/{remote_rel_path}"
                    files_to_delete_list.append(full_remote_path)

            # Delete files in parallel if any exist
            if files_to_delete_list:
                delete_success, delete_results = self.delete_files(
                    files_to_delete_list, progress_callback=progress_callback
                )
                if not delete_success:
                    success = False
                # Merge deletion results into the main status dictionary
                file_statuses.update(delete_results)

        return success, file_statuses

    def check_health(self) -> bool:
        """
        Check if the file service is healthy.
        This method uses the health endpoint which doesn't require authentication.
        Returns:
            True if the service is healthy, False otherwise
        """
        try:
            response = self.client.httpx_request_and_process(genapi._get_check_kwargs())

            # Check if the response indicates the service is healthy
            return response.status_code == 200 and response.parsed_json is not None
        except Exception as e:
            self.logger.error(f"Health check failed: {e}")
            return False

    def check_auth_ok(self) -> bool:
        """
        Check if authentication is working correctly.
        This method attempts to list objects, which requires authentication.
        Returns:
            True if authentication is working, False otherwise
        """
        try:
            # Try to list objects in a default location

            call_args = genapi._get_get_object_kwargs(
                orgFriendlyId=self.organization_id,
                filePath="",
                list_type=2,  # Default list type
            )
            response = self.client.httpx_request_and_process(call_args)

            # If we get a 200 or 404, authentication is working
            # 404 would mean the path doesn't exist, but auth is still valid
            return response.status_code in (200, 404)
        except Exception as e:
            self.logger.error(f"Authentication check failed: {e}")
            return False


def _is_notebook() -> bool:
    """Determine if code is running in a Jupyter notebook environment."""
    try:
        shell = get_ipython().__class__.__name__
        # If shell is either ZMQInteractiveShell (Jupyter notebook or qtconsole)
        # or TerminalInteractiveShell (IPython terminal), we return True/False
        if shell == "ZMQInteractiveShell":
            return True  # Jupyter notebook or qtconsole
        elif shell == "TerminalInteractiveShell":
            return False
    except NameError:
        pass  # Standard Python interpreter
    return False


def create_files_tqdm_callback(description: str) -> Callable[[int, int], None] | None:
    """
    Creates a callback closure to manage a tqdm progress bar, usable to pass to
    download_files, sync_folder_up, etc. First callback initializes total.
    """
    pbar = None
    try:
        # Import appropriate tqdm version based for runing in notebook/not.
        if _is_notebook():
            from tqdm.notebook import tqdm
        else:
            from tqdm.auto import tqdm

        def _callback(completed, total, status=None):
            nonlocal pbar
            if total <= 0:
                return
            if pbar is None:
                pbar = tqdm(total=total, desc=description)

            # Update progress count
            pbar.n = completed
            # If status information is provided, update the postfix
            if status:
                if isinstance(status, str):
                    pbar.set_postfix_str(status)
                elif isinstance(status, dict):
                    pbar.set_postfix(**status)

            pbar.refresh()
            if completed >= total:
                pbar.close()
                pbar = None

        return _callback
    except ImportError:
        return None
