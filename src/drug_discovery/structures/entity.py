"""
This module defines the Entity class for handling file uploads to a remote server in the context of drug discovery structures.

The Entity class provides methods to manage and upload files, such as protein structure files, to a remote storage system using the DeepOrigin FilesClient.
"""

from dataclasses import dataclass, field
import os
from typing import Optional

from deeporigin.files import FilesClient


@dataclass
class Entity:
    """
    Represents an entity with file upload capabilities to a remote server.

    This class manages the remote path and provides an upload method to ensure that the entity's file is uploaded to the remote storage if it does not already exist there. It uses the DeepOrigin FilesClient for remote file operations.
    """

    _remote_path: Optional[str] = field(default=None, init=False)

    def upload(self):
        """Upload the entity to the remote server."""

        files_client = FilesClient()

        remote_files = files_client.list_folder(self._remote_path_base, recursive=True)
        remote_files = list(remote_files.keys())

        files_to_upload = {}

        protein_path = self._remote_path_base + os.path.basename(self.file_path)
        if protein_path not in remote_files:
            files_to_upload[str(self.file_path)] = protein_path

        # set protein path
        self._remote_path = protein_path
