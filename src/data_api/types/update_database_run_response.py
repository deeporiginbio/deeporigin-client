# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.



from .._models import BaseModel
from .database import Database

__all__ = ["UpdateDatabaseRunResponse"]


class UpdateDatabaseRunResponse(BaseModel):
    data: Database
