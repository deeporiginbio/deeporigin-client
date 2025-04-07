from enum import Enum


class SyncFileSchemaDtoProvider(str, Enum):
    S3 = "S3"

    def __str__(self) -> str:
        return str(self.value)
