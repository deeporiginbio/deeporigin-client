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
    files = client.list_dir('files:///data/', flag="recursive")
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

import logging
import os
import re
from typing import Any, Dict, List, Optional, Union

from deeporigin import auth
from deeporigin.utils.network import _get_domain_name

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

# from .file_service.models.sync_file_schema_dto import SyncFileSchemaDto
# from .file_service.models.sync_file_schema_dto_credentials import (
#    SyncFileSchemaDtoCredentials,
# )
# from .file_service.models.sync_file_schema_dto_provider import SyncFileSchemaDtoProvider


logger = logging.getLogger(__name__)


class FilesClient:
    """
    Client for interacting with the DeepOrigin file service.

    This client provides high-level methods for file operations like upload, download,
    listing directories, and synchronizing files between local and remote storage.
    """

    FILES_URL_PATTERN = re.compile(r"^files:///(?P<path>.*)$")

    def __init__(
        self,
        base_url: Optional[str] = None,
        token: Optional[str] = None,
        verify_ssl: Union[bool, str] = True,
    ):
        """
        Initialize the FilesClient.
        Args:
            base_url: Base URL of the file service API
            token: Authentication token (optional)
            verify_ssl: Whether to verify SSL certificates or path to CA bundle
        """

        if not base_url:
            base_url = _get_domain_name()
        self._base_url = base_url
        if not token:
            tokens = auth.get_tokens(refresh=False)
            token = tokens["access"]

        if token:
            self.client = AuthenticatedClient(
                base_url=self._base_url, token=token, verify_ssl=verify_ssl
            )
        else:
            self.client = Client(base_url=self._base_url, verify_ssl=verify_ssl)

    def get_base_url(self) -> str:
        """
        Get the base URL of the file service API, mostly for visibilty purposes.
        Returns:
            Base URL of the file service API
        """
        return self._base_url

    def _parse_files_url(self, url: str) -> tuple[str, str]:
        """
        Parse a files:/// URL to extract org ID and path.
        Args:
            url: URL in the format files:///path
        Returns:
            Tuple of (org_friendly_id, path)
        Raises:
            ValueError: If the URL is not in the expected format
        """
        match = self.FILES_URL_PATTERN.match(url)
        if not match:
            raise ValueError(
                f"Invalid files URL: {url}. Expected format: files:///path"
            )

        path = match.group("path")

        # In a real implementation, we might extract the org_friendly_id from the path
        # For now, we'll use a default value
        org_friendly_id = "default"

        return org_friendly_id, path

    def list_dir(self, path: str, flag: Optional[str] = None) -> List[Dict[str, Any]]:
        """
        List files and directories at the specified path.
        Args:
            path: Path in the format files:///path
            flag: Optional flags, e.g., "recursive"
        Returns:
            List of file/directory metadata objects
        """
        org_friendly_id, remote_path = self._parse_files_url(path)

        # The list_type parameter seems to be required in the API
        # We'll map the flag to an appropriate value
        list_type = 2.0  # Default value
        if flag == "recursive":
            list_type = 1.0

        response = get_object.sync_detailed(
            org_friendly_id=org_friendly_id,
            file_path=remote_path,
            client=self.client,
            list_type=list_type,
        )

        #  print("response: ", response)

        if response.status_code == 200:
            # Parse the response content as needed
            # This will depend on the actual response format from the API
            return self._parse_list_response(response.content)
        else:
            logger.error(f"Failed to list directory {path}: {response.status_code}")
            response.raise_for_status()
            return []

    def _parse_list_response(self, content: bytes) -> List[Dict[str, Any]]:
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
        """
        Check if a file exists in remote storage.

        Args:
            path: Path in the format files:///path

        Returns:
            True if the file exists, False otherwise
        """
        try:
            metadata = self.get_metadata(path)
            # If we got metadata without an error, the file exists
            return bool(metadata)
        except Exception:
            # If there was an error getting metadata, the file probably doesn't exist
            return False

    def upload_file(self, src: str, dest: str, overwrite: bool = False) -> bool:
        """
        Upload a file from local path to remote storage.

        Args:
            src: Local file path
            dest: Remote path in the format files:///path
            overwrite: If True, delete the file first if it exists

        Returns:
            True if successful, False otherwise
        """
        if not os.path.exists(src):
            raise FileNotFoundError(f"Source file not found: {src}")

        # Check if file exists and delete if overwrite is True
        if overwrite and self._remote_file_exists(dest):
            logger.info(f"Deleting existing file at {dest} before upload")
            if not self.delete_file(dest):
                logger.warning(f"Failed to delete existing file at {dest}")
                # Continue with upload anyway

        org_friendly_id, remote_path = self._parse_files_url(dest)

        # Create a File object first
        with open(src, "rb") as f:
            file_obj = File(
                payload=f,
                file_name=os.path.basename(src),
                mime_type="application/octet-stream",
            )

            # Create the PutObjectBody with the File object
            body = PutObjectBody(file=file_obj)

            # Send the request with the proper body and include the remote path
            response = put_object.sync_detailed(
                org_friendly_id=org_friendly_id,
                file_path=remote_path,
                client=self.client,
                body=body,
            )

        return response.status_code == 200

    def download_file(self, src: str, dest: str) -> bool:
        """
        Download a file from remote storage to local path.
        Args:
            src: Remote path in the format files:///path
            dest: Local file path
        Returns:
            True if successful, False otherwise
        """
        org_friendly_id, remote_path = self._parse_files_url(src)

        # The get_object API will need to be called differently for file download
        # vs. directory listing. This implementation assumes that omitting list_type
        # will result in a file download.
        response = get_object.sync_detailed(
            org_friendly_id=org_friendly_id, file_path=remote_path, client=self.client
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

    def get_metadata(self, path: str) -> Dict[str, Any]:
        """
        Get metadata about a file or directory.
        Args:
            path: Path in the format files:///path
        Returns:
            Metadata dictionary
        """
        org_friendly_id, remote_path = self._parse_files_url(path)

        response = head_object.sync_detailed(
            org_friendly_id=org_friendly_id, file_path=remote_path, client=self.client
        )

        if response.status_code == 200:
            # Extract metadata from response headers
            return dict(response.headers)
        else:
            logger.error(f"Failed to get metadata for {path}: {response.status_code}")
            return {}

    def sync_dir(self, src: str, dst: str) -> bool:
        """
        Synchronize files between local and remote directories.
        Args:
            src: Source path (can be local or files:/// URL)
            dst: Destination path (can be local or files:/// URL)
        Returns:
            True if successful, False otherwise
        """
        # Determine sync direction
        if src.startswith("files:///") and not dst.startswith("files:///"):
            # Download from remote to local
            return self._sync_remote_to_local(src, dst)
        elif not src.startswith("files:///") and dst.startswith("files:///"):
            # Upload from local to remote
            return self._sync_local_to_remote(src, dst)
        else:
            raise ValueError("Either src or dst must be a files:/// URL, but not both")

    def _sync_remote_to_local(self, remote_path: str, local_path: str) -> bool:
        """
        Sync files from remote to local directory.
        Args:
            remote_path: Remote path in the format files:///path
            local_path: Local directory path
        Returns:
            True if successful, False otherwise
        """
        # Get list of remote files with metadata
        remote_files = self.list_dir(remote_path, flag="recursive")

        # Ensure local directory exists
        os.makedirs(local_path, exist_ok=True)

        success = True
        for file_info in remote_files:
            # Extract the relative path from the Key
            rel_path = file_info.get("Key", "")

            # Don't append the Key to remote_path - it already contains the full path
            # Instead, construct a new files:/// URL with the path
            remote_file_path = f"files:///{rel_path}"

            # For the local path, join the base local path with just the filename part
            local_file_path = os.path.join(local_path, os.path.basename(rel_path))

            # If the file is in a subdirectory, make sure that directory exists locally
            local_dir = os.path.dirname(local_file_path)
            if local_dir:
                os.makedirs(local_dir, exist_ok=True)

            # Check if we need to download this file
            if os.path.exists(local_file_path):
                local_mtime = os.path.getmtime(local_file_path)
                # Convert the remote timestamp to a comparable format
                remote_mtime_str = file_info.get("LastModified", "")
                try:
                    # Parse ISO format timestamp
                    from datetime import datetime

                    remote_mtime = datetime.fromisoformat(
                        remote_mtime_str.replace("Z", "+00:00")
                    ).timestamp()
                except (ValueError, TypeError):
                    remote_mtime = 0

                if local_mtime >= remote_mtime:
                    continue  # Local file is newer or same age

            # Download the file
            if not self.download_file(remote_file_path, local_file_path):
                logger.error(
                    f"Failed to download file {remote_file_path} to {local_file_path}"
                )
                success = False

        return success

    def _sync_local_to_remote(self, local_path: str, remote_path: str) -> bool:
        """
        Sync files from local to remote directory.
        Args:
            local_path: Local directory path
            remote_path: Remote path in the format files:///path
        Returns:
            True if successful, False otherwise
        """
        # Get list of remote files with metadata to compare
        try:
            remote_files = {
                file_info["Key"]: file_info.get("LastModified", 0)
                for file_info in self.list_dir(remote_path, flag="recursive")
            }
        except Exception as e:
            logger.warning(f"Could not list remote directory, assuming it's empty: {e}")
            remote_files = {}

        success = True
        for root, _, files in os.walk(local_path):
            for filename in files:
                local_file_path = os.path.join(root, filename)

                # Get relative path from the local_path
                rel_path = os.path.relpath(local_file_path, local_path)
                remote_file_path = f"{remote_path}/{rel_path}"

                # Check if we need to upload this file
                local_mtime = os.path.getmtime(local_file_path)
                remote_mtime = remote_files.get(rel_path, 0)

                if local_mtime <= remote_mtime:
                    continue  # Remote file is newer or same age

                # Upload the file
                if not self.upload_file(local_file_path, remote_file_path):
                    success = False

        return success

    def delete_file(self, path: str) -> bool:
        """
        Delete a file from remote storage.
        Args:
            path: Path in the format files:///path
        Returns:
            True if successful, False otherwise
        """
        org_friendly_id, remote_path = self._parse_files_url(path)

        response = delete_object.sync_detailed(
            org_friendly_id=org_friendly_id, file_path=remote_path, client=self.client
        )

        return response.status_code == 200

    '''
    def sync_with_external_provider(self, path: str, provider_config: Dict[str, Any]) -> bool:
        """
        Sync files with an external provider (e.g., S3).        
        Args:
            path: Path in the format files:///path
            provider_config: Provider configuration
            
        Returns:
            True if successful, False otherwise
        """
        org_friendly_id, remote_path = self._parse_files_url(path)
        
        # Create the sync request body
        credentials = SyncFileSchemaDtoCredentials(
            region=provider_config.get("region", ""),
            secret_key=provider_config.get("secretKey", ""),
            access_key=provider_config.get("accessKey", "")
        )
        
        sync_request = SyncFileSchemaDto(
            provider=SyncFileSchemaDtoProvider.S3,
            credentials=credentials,
            path=provider_config.get("path", "")
        )
        
        response = sync_objects.sync_detailed(
            org_friendly_id=org_friendly_id,
            file_path=remote_path,
            client=self.client,
            json_body=sync_request
        )
        
        return response.status_code == 201
    '''

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
                list_type=2.0,  # Default list type
            )

            # If we get a 200 or 404, authentication is working
            # 404 would mean the path doesn't exist, but auth is still valid
            return response.status_code in (200, 404)
        except Exception as e:
            logger.error(f"Authentication check failed: {e}")
            return False
