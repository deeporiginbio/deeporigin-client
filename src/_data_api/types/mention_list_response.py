# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from typing import List, Optional
from typing_extensions import Literal

from pydantic import Field as FieldInfo

from .._models import BaseModel

__all__ = ["MentionListResponse", "Data", "DataMention"]


class DataMention(BaseModel):
    id: str

    hid: str

    type: Literal["database", "row", "workspace"]

    name: Optional[str] = None

    parent_id: Optional[str] = FieldInfo(alias="parentId", default=None)


class Data(BaseModel):
    mentions: List[DataMention]


class MentionListResponse(BaseModel):
    data: Data
