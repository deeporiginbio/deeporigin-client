# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from typing import List, Optional
from typing_extensions import Literal

from pydantic import Field as FieldInfo

from .._models import BaseModel

__all__ = ["Database", "Col"]


class Col:
    pass


class Database(BaseModel):
    id: str

    date_created: str = FieldInfo(alias="dateCreated")

    hid: str

    hid_prefix: Optional[str] = FieldInfo(alias="hidPrefix", default=None)

    name: Optional[str] = None

    type: Literal["database"]

    cols: Optional[List[Col]] = None

    created_by_user_drn: Optional[str] = FieldInfo(alias="createdByUserDrn", default=None)

    creation_block_id: Optional[str] = FieldInfo(alias="creationBlockId", default=None)

    creation_parent_id: Optional[str] = FieldInfo(alias="creationParentId", default=None)

    date_updated: Optional[str] = FieldInfo(alias="dateUpdated", default=None)

    edited_by_user_drn: Optional[str] = FieldInfo(alias="editedByUserDrn", default=None)

    editor: Optional[object] = None

    is_template: Optional[bool] = FieldInfo(alias="isTemplate", default=None)

    parent_id: Optional[str] = FieldInfo(alias="parentId", default=None)

    row_json_schema: Optional[object] = FieldInfo(alias="rowJsonSchema", default=None)

    submission_status: Optional[Literal["draft", "final"]] = FieldInfo(alias="submissionStatus", default=None)

    validation_status: Optional[str] = FieldInfo(alias="validationStatus", default=None)
