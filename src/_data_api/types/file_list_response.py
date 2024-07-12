# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from typing import List, Optional

from pydantic import Field as FieldInfo

from .file import File
from .._models import BaseModel

__all__ = ["FileListResponse", "Data", "DataAssignment"]


class DataAssignment(BaseModel):
    row_id: str = FieldInfo(alias="rowId")


class Data(BaseModel):
    file: File

    assignments: Optional[List[DataAssignment]] = None


class FileListResponse(BaseModel):
    data: List[Data]
