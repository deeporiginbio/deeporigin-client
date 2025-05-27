from dataclasses import dataclass, field
import os
from typing import Optional

from deeporigin.files import FilesClient


@dataclass
class Entity:
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
