"""
This module defines the Entity class for handling file uploads to a remote server in the context of drug discovery structures.

The Entity class provides methods to manage and upload files, such as protein structure files, to a remote storage system using the DeepOrigin FilesClient.
"""

from dataclasses import dataclass, field
import os
from typing import Optional


@dataclass
class Entity:
    """
    Represents an entity with file upload capabilities to a remote server.

    This class manages the remote path and provides an upload method to ensure that the entity's file is uploaded to the remote storage if it does not already exist there. It uses the DeepOrigin FilesClient for remote file operations.
    """

    _remote_path: Optional[str] = field(default=None, init=False)

    def upload(self):
        """Upload the entity to the remote server."""

        from deeporigin.platform import file_api

        remote_files = file_api.get_object_directory(
            file_path="/entities/",
            recursive=True,
        )
        remote_files = [file.Key for file in remote_files]

        remote_path = self._remote_path_base + os.path.basename(self.file_path)
        if remote_path not in remote_files:
            file_api.upload_file(remote_path=remote_path, local_path=self.file_path)

        # set protein path
        self._remote_path = remote_path
