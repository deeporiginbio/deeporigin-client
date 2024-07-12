# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.



from .._models import BaseModel
from .database import Database

__all__ = ["UpdateDatabaseColumnRunResponse", "Data", "DataColumn"]


class DataColumn:
    pass


class Data(BaseModel):
    column: DataColumn

    database: Database


class UpdateDatabaseColumnRunResponse(BaseModel):
    data: Data
