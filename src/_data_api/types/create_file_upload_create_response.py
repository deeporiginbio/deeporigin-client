# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.


from pydantic import Field as FieldInfo

from .file import File
from .._models import BaseModel

__all__ = ["CreateFileUploadCreateResponse", "Data"]


class Data(BaseModel):
    file: File

    upload_url: str = FieldInfo(alias="uploadUrl")


class CreateFileUploadCreateResponse(BaseModel):
    data: Data
