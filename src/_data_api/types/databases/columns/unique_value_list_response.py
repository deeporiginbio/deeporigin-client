# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from typing import List

from ...._models import BaseModel

__all__ = ["UniqueValueListResponse", "Data"]


class Data(BaseModel):
    name: str

    value: str


class UniqueValueListResponse(BaseModel):
    data: List[Data]
