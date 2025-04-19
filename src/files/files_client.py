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
        src="files:///path/to/remote/file.txt",
        dest="/local/path/file.txt"
    )

    if success:
        print("File downloaded successfully")
    else:
        print("Download failed")
    ```
"""

from dataclasses import dataclass, field
from datetime import datetime
from enum import Enum
import logging
import os
import re
import concurrent.futures
import time

from deeporigin import auth
from deeporigin.utils.network import _get_domain_name
from deeporigin.utils.config import get_value

from .file_service import AuthenticatedClient, Client
from .file_service.api.default import (
    delete_object,
    get_object,
    head_object,
    put_object,
    # sync_objects,
)
from .file_service.models.put_object_body import PutObjectBody
from .file_service.types import File

logger = logging.getLogger(__name__)

# Remote path prefix for file service URLs
#REMOTE_PATH_PREFIX = "files:///"
#FILES_URL_PATTERN = re.compile(rf"^{REMOTE_PATH_PREFIX}(?P<path>.*)$")
REMOTE_PATH_PREFIX = ""


class DirSyncMode(Enum):
    """
    Enum for directory synchronization modes.
    """

    # We don't know if the destination file not present locally was deleted locally or not.
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
    def default(cls) -> "DirSyncMode":
        """Return the default sync mode (REPLACE)"""
        return cls.REPLACE

    def can_delete_nonexisting(self) -> bool:
        return self == DirSyncMode.REPLACE or self == DirSyncMode.FORCE_OVERWRITE


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
    listing directories, and synchronizing files between local and remote storage.
    """

    # FilesClient instance attributes
    base_url: str
    organization_id: str
    client: AuthenticatedClient | Client
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

        self.io_max_retries = 3
        self.io_retry_delay = 0.5

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
        self.client = AuthenticatedClient(
            base_url=self.base_url, token=token, verify_ssl=verify_ssl
        )


    def _extract_path_from_url(self, url: str) -> str:
        # This was used in the past to extract the path from the URL, seems unnecessary.
        # match = FILES_URL_PATTERN.match(url)
        # return match.group("path")
        return url

    def list_folder(self, path: str, flag: str | None = None) -> dict[str, FileMetadata]:
        """
        List files and folders at the specified path. Returned dictionary keys will have
        the full remote path as the key, including the provided path argument part.
        This means that the resulting path key can be used directly for download, etc.
        Args:
            path: Path in the format files:///path
            flag: Optional flag for listing behavior
        Returns:
            Dictionary mapping file paths to FileMetadata objects, maintaining the order
            from the response
        """
        remote_path = self._extract_path_from_url(path)

        if flag == "recursive":
            raise ValueError("list_folder recursive flag not implemented")

        # The list_type parameter seems to be required in the API
        # We'll map the flag to an appropriate value
        list_type = 2  # Default value

        response = get_object.sync_detailed(
            org_friendly_id=self.organization_id,
            file_path=remote_path,
            client=self.client,
            list_type=list_type,
        )

        #  print("response: ", response)

        if response.status_code == 200:
            # Parse the response content as needed
            # This will depend on the actual response format from the API
            raw_files = self._parse_list_response(response.content)
            
            # Create a dictionary with key_path as keys, maintaining order
            result = {}
            for file_data in raw_files:
                metadata = FileMetadata.from_dict(
                    file_data, mapping=FileMetadata.METADATA_LIST_FIELD_MAPPING
                )
                if metadata.key_path:
                    result[metadata.key_path] = metadata
            
            return result
        else:
            logger.error(f"Failed to list directory {path}: {response.status_code}")
            response.raise_for_status()
            return {}

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
            logger.warning("Could not parse list response as JSON")
            return []

    def _remote_file_exists(self, path: str) -> bool:
        # If we got metadata without an error, the file exists
        try:
            metadata = self.get_metadata(path)            
            return bool(metadata)
        except Exception:
            # If there was an error getting metadata, the file probably doesn't exist
            return False

    def upload_file(self, src: str, dest: str) -> bool:
        """
        Upload a file from local path to remote storage.
        Args:
            src: Local source path
            dest: Remote destination path
            overwrite: Whether to overwrite existing file (default: False)
        Returns:
            True if upload was successful, False otherwise
        """
        # Use the new upload_files method with a single file
        success, _ = self.upload_files({src: dest})
        return success

    def upload_files(self, src_to_dest: dict[str, str]) -> tuple[bool, dict[str, str]]:
        """
        Upload multiple files from local paths to remote storage in parallel.
        Args:
            src_to_dest: Dictionary mapping source file paths to destination paths
        Returns:
            Tuple containing:
                - Boolean indicating if all uploads were successful
                - Dictionary mapping source paths to status ("OK" or error message)
        """
        if not src_to_dest:
            return True, {}

        results = {}
        all_success = True

        # For a single file, use _upload_single_file directly
        if len(src_to_dest) == 1:
            src, dest = next(iter(src_to_dest.items()))
            success, status = self._upload_single_file(src, dest)
            results[src] = status
            return success, results
        else:
            # For multiple files, use concurrent futures
            with concurrent.futures.ThreadPoolExecutor() as executor:
                future_to_src = {}
                
                for src, dest in src_to_dest.items():
                    if not os.path.exists(src):
                        results[src] = f"Source file not found: {src}"
                        all_success = False
                        continue
                        
                    future = executor.submit(self._upload_single_file, src, dest)
                    future_to_src[future] = src
                
                for future in concurrent.futures.as_completed(future_to_src):
                    src = future_to_src[future]
                    try:
                        success, status = future.result()
                        results[src] = status
                        if not success:
                            all_success = False
                    except Exception as e:
                        results[src] = f"Upload failed with error: {str(e)}"
                        all_success = False

        return all_success, results

    def _upload_single_file(self, src: str, dest: str) -> tuple[bool, str]:
        remote_path = self._extract_path_from_url(dest)
        last_error = None
        
        for attempt in range(self.io_max_retries):
            try:
                if not os.path.exists(src):                    
                    return False,  f"Source file not found: {src}"

                with open(src, "rb") as f:
                    file_obj = File(
                        payload=f,
                        file_name=os.path.basename(src),
                        mime_type="application/octet-stream",
                    )

                    body = PutObjectBody(file=file_obj)

                    response = put_object.sync_detailed(
                        org_friendly_id=self.organization_id,
                        file_path=remote_path,
                        client=self.client,
                        body=body,
                    )

                if response.status_code == 200:
                    return True, "OK"
                
                # Determine if we should retry based on status code. 500 codes are server errors.
                # 429 is rate limit exceeded, 408 is request timeout.
                if response.status_code in [429, 408] or (500 <= response.status_code < 600):
                    if attempt < self.io_max_retries - 1:
                        last_error = f"Attempt {attempt+1} failed with status {response.status_code}"
                        logger.warning(f"{last_error}, retrying {src} upload in {self.io_retry_delay}s...")
                        time.sleep(self.io_retry_delay)
                        continue
                
                # For other status codes, don't retry
                return False, f"Upload failed with status code: {response.status_code}"
                
            except (ConnectionError, TimeoutError) as e:
                # Retry on connection errors
                if attempt < self.io_max_retries - 1:
                    last_error = f"Connection error on attempt {attempt+1}: {str(e)}"
                    logger.warning(f"{last_error}, retrying {src} upload in {self.io_retry_delay}s")
                    time.sleep(self.io_retry_delay)
                    continue
                return False, f"Upload failed with connection error: {str(e)}"
                
            except Exception as e:
                return False, f"Upload failed with error: {str(e)}"
        
        # We've exhausted all retries, return the last error message
        return False, last_error


    def download_file(self, src: str, dest: str) -> bool:
        """
        Download a file from remote storage to local path.
        Args:
            src: Remote source path in the format files:///path
            dest: Local destination path
        Returns:
            True if download was successful, False otherwise
        """
        remote_path = self._extract_path_from_url(src)

        # The get_object API will need to be called differently for file download
        # vs. directory listing. This implementation assumes that omitting list_type
        # will result in a file download.
        response = get_object.sync_detailed(
            org_friendly_id=self.organization_id,
            file_path=remote_path,
            client=self.client,
        )  # list_type=0.0  - wasn't working # Assuming 0.0 means download file
        # looks like bad if(type) check in file-service controller --TBD Niels

        if response.status_code == 200:
            # Ensure the directory exists
            os.makedirs(os.path.dirname(os.path.abspath(dest)), exist_ok=True)

            # Write the content to the destination file
            with open(dest, "wb") as f:
                f.write(response.content)
            return True
        else:
            logger.error(f"Failed to download file {src}: {response.status_code}")
            return False

    def get_metadata(self, path: str) -> FileMetadata | None:
        """
        Get metadata about a file or directory.
        Args:
            path: Path in the format files:///path
        Returns:
            FileMetadata object containing the metadata
        """
        remote_path = self._extract_path_from_url(path)

        response = head_object.sync_detailed(
            org_friendly_id=self.organization_id,
            file_path=remote_path,
            client=self.client,
        )

        if response.status_code == 200:
            # Extract metadata from response headers and create a FileMetadata object
            return FileMetadata.from_dict(dict(response.headers), key_path=remote_path)
        else:
            logger.error(f"Failed to get metadata for {path}: {response.status_code}")
            return None


    def sync_folder_down(
        self, remote_path: str, local_path: str, mode: DirSyncMode = DirSyncMode.REPLACE
    ) -> tuple[bool, dict[str, str]]:
        """
        Sync files from remote to local directory.
        Args:
            remote_path: Remote path in the format files:///path
            local_path: Local directory path
            mode: Synchronization mode
        Returns:
            Tuple containing:
                - Boolean indicating if the sync was successful
                - Dictionary mapping file paths to their status
        """
        # Get list of remote files with metadata
        remote_files_dict = self.list_folder(remote_path)  # todo: add recursive flag

        # Extract the base remote path without the prefix for path comparisons
        base_remote_path = self._extract_path_from_url(remote_path)

        # Ensure local directory exists
        os.makedirs(local_path, exist_ok=True)

        # Track all file operations
        file_statuses = {}
        success = True

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
                if mode == DirSyncMode.FORCE_OVERWRITE:
                    should_update = True
                elif mode == DirSyncMode.REPLACE:
                    should_update = local_mtime < remote_mtime
                elif mode == DirSyncMode.MERGE:
                    should_update = local_mtime < remote_mtime
                elif mode == DirSyncMode.COPY_NON_EXISTING:
                    should_update = False

                if should_update:
                    # Download the file
                    if self.download_file(full_remote_file_path, local_file_path):
                        file_statuses[local_file_path] = "copied"
                    else:
                        logger.error(
                            f"Failed to download file {full_remote_file_path} to {local_file_path}"
                        )
                        file_statuses[local_file_path] = "error"
                        success = False
                else:
                    file_statuses[local_file_path] = "kept"
            else:
                # File doesn't exist locally, download it
                if self.download_file(full_remote_file_path, local_file_path):
                    file_statuses[local_file_path] = "copied"
                else:
                    logger.error(
                        f"Failed to download file {full_remote_file_path} to {local_file_path}"
                    )
                    file_statuses[local_file_path] = "error"
                    success = False

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
                        file_statuses[full_local_path] = "deleted"
                    except Exception as e:
                        logger.error(f"Failed to delete file {full_local_path}: {e}")
                        file_statuses[full_local_path] = "error"
                        success = False

        return success, file_statuses

    def sync_folder_up(
        self, local_path: str, remote_path: str, mode: DirSyncMode = DirSyncMode.REPLACE
    ) -> tuple[bool, dict[str, str]]:
        """
        Sync files from local to remote directory.
        Args:
            local_path: Local directory path
            remote_path: Remote path in the format files:///path
            mode: Synchronization mode
        Returns:
            Tuple containing:
                - Boolean indicating if the sync was successful
                - Dictionary mapping file paths to their status
        """
        # Extract the base remote path without the prefix for path comparisons
        base_remote_path = self._extract_path_from_url(remote_path)

        # Get list of remote files with metadata to compare
        remote_files_dict = {}
        try:
            remote_files_dict = self.list_folder(remote_path)  # todo: add recursive flag
        except Exception as e:
            logger.warning(f"Could not list remote directory, assuming it's empty: {e}")            

        # Track all file operations
        file_statuses = {}
        success = True

        # Process local files
        for root, _, files in os.walk(local_path):
            for filename in files:
                local_file_path = os.path.join(root, filename)

                # Get relative path from the local_path
                rel_path = os.path.relpath(local_file_path, local_path)
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
                    if mode == DirSyncMode.FORCE_OVERWRITE:
                        should_update = True
                    elif mode == DirSyncMode.REPLACE:
                        should_update = local_mtime > remote_mtime
                    elif mode == DirSyncMode.MERGE:
                        should_update = local_mtime > remote_mtime
                    elif mode == DirSyncMode.COPY_NON_EXISTING:
                        should_update = False

                    if should_update:
                        # Upload the file
                        if self.upload_file(local_file_path, remote_file_path):
                            file_statuses[remote_file_path] = "copied"
                        else:
                            logger.error(
                                f"Failed to upload file {local_file_path} to {remote_file_path}"
                            )
                            file_statuses[remote_file_path] = "error"
                            success = False
                    else:
                        file_statuses[remote_file_path] = "kept"
                else:
                    # File doesn't exist remotely, upload it
                    if self.upload_file(local_file_path, remote_file_path):
                        file_statuses[remote_file_path] = "copied"
                    else:
                        logger.error(
                            f"Failed to upload file {local_file_path} to {remote_file_path}"
                        )
                        file_statuses[remote_file_path] = "error"
                        success = False

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
            for remote_file_key in remote_files_dict:
                if (
                    os.path.relpath(remote_file_key, base_remote_path)
                    not in local_files
                ):
                    full_remote_path = f"{remote_path}/{remote_file_key}"
                    if self.delete_file(full_remote_path):
                        file_statuses[full_remote_path] = "deleted"
                    else:
                        logger.error(f"Failed to delete file {full_remote_path}")
                        file_statuses[full_remote_path] = "error"
                        success = False

        return success, file_statuses

    def delete_file(self, path: str) -> bool:
        """
        Delete a file from remote storage.
        Args:
            path: Path in the format files:///path
        Returns:
            True if deletion was successful, False otherwise
        """
        remote_path = self._extract_path_from_url(path)

        response = delete_object.sync_detailed(
            org_friendly_id=self.organization_id,
            file_path=remote_path,
            client=self.client,
        )

        return response.status_code == 200

    def check_health(self) -> bool:
        """
        Check if the file service is healthy.
        This method uses the health endpoint which doesn't require authentication.
        Returns:
            True if the service is healthy, False otherwise
        """
        try:
            # Create a non-authenticated client with the same parameters
            non_auth_client = Client(
                base_url=self.client._base_url, verify_ssl=self.client._verify_ssl
            )

            # Call the health check endpoint
            from .file_service.api.health import check

            response = check.sync_detailed(client=non_auth_client)

            # Check if the response indicates the service is healthy
            return response.status_code == 200 and response.parsed is not None
        except Exception as e:
            logger.error(f"Health check failed: {e}")
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
            org_friendly_id = "default"

            # Use get_object with list_type to check authentication
            response = get_object.sync_detailed(
                org_friendly_id=org_friendly_id,
                file_path="",
                client=self.client,
                list_type=2,  # Default list type
            )

            # If we get a 200 or 404, authentication is working
            # 404 would mean the path doesn't exist, but auth is still valid
            return response.status_code in (200, 404)
        except Exception as e:
            logger.error(f"Authentication check failed: {e}")
            return False
