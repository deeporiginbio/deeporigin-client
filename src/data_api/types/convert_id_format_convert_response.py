# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from typing import List

from .._models import BaseModel

__all__ = ["ConvertIDFormatConvertResponse", "Data"]


class Data(BaseModel):
    id: str

    hid: str


class ConvertIDFormatConvertResponse(BaseModel):
    data: List[Data]
