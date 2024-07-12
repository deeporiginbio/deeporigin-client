# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from typing import List, Optional
from typing_extensions import Literal

from .._models import BaseModel

__all__ = ["SeqData", "Annotation"]


class Annotation(BaseModel):
    end: float

    name: str

    start: float

    color: Optional[str] = None

    direction: Optional[Literal[1, 0, -1]] = None

    type: Optional[str] = None


class SeqData(BaseModel):
    seq: str

    annotations: Optional[List[Annotation]] = None

    name: Optional[str] = None

    type: Optional[Literal["dna", "rna", "aa", "unknown"]] = None
