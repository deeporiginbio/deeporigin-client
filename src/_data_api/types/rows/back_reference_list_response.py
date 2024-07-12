# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from typing import List

from ..._models import BaseModel

__all__ = ["BackReferenceListResponse", "Data", "DataRow"]


class DataRow:
    pass


class Data(BaseModel):
    rows: List[DataRow]


class BackReferenceListResponse(BaseModel):
    data: Data
