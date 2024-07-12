# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from typing import Optional

from pydantic import Field as FieldInfo

from .._models import BaseModel

__all__ = ["File"]


class File(BaseModel):
    id: str

    content_length: float = FieldInfo(alias="contentLength")

    content_type: Optional[str] = FieldInfo(alias="contentType", default=None)

    created_by_user_drn: Optional[str] = FieldInfo(
        alias="createdByUserDrn", default=None
    )

    date_created: str = FieldInfo(alias="dateCreated")

    date_updated: str = FieldInfo(alias="dateUpdated")

    name: str

    status: str

    uri: str
