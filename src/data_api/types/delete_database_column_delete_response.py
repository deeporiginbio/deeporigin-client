# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.



from .._models import BaseModel
from .database import Database

__all__ = ["DeleteDatabaseColumnDeleteResponse", "Data"]


class Data(BaseModel):
    database: Database


class DeleteDatabaseColumnDeleteResponse(BaseModel):
    data: Data
