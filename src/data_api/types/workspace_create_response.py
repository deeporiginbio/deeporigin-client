# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.


from .._models import BaseModel
from .workspace import Workspace

__all__ = ["WorkspaceCreateResponse"]


class WorkspaceCreateResponse(BaseModel):
    data: Workspace
