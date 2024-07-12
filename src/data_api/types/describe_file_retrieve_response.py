# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.


from .file import File
from .._models import BaseModel

__all__ = ["DescribeFileRetrieveResponse"]


class DescribeFileRetrieveResponse(BaseModel):
    data: File
