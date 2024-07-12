# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from typing import List, Optional

from pydantic import Field as FieldInfo

from .._models import BaseModel

__all__ = ["ConfigureColumnSelectOptionConfigureResponse", "Data", "DataConfigSelect"]


class DataConfigSelect(BaseModel):
    options: List[str]

    can_create: Optional[bool] = FieldInfo(alias="canCreate", default=None)


class Data(BaseModel):
    id: str

    config_select: DataConfigSelect = FieldInfo(alias="configSelect")


class ConfigureColumnSelectOptionConfigureResponse(BaseModel):
    data: Data
