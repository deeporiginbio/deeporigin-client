# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing import Union, Iterable
from typing_extensions import Required, TypedDict

__all__ = ["ConvertIDFormatConvertParams", "Conversion", "ConversionID", "ConversionHid"]


class ConvertIDFormatConvertParams(TypedDict, total=False):
    conversions: Required[Iterable[Conversion]]


class ConversionID(TypedDict, total=False):
    id: Required[str]


class ConversionHid(TypedDict, total=False):
    hid: Required[str]


Conversion = Union[ConversionID, ConversionHid]
