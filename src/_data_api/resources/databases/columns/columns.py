# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from ...._compat import cached_property
from ...._resource import SyncAPIResource, AsyncAPIResource
from .unique_values import (
    UniqueValuesResource,
    AsyncUniqueValuesResource,
    UniqueValuesResourceWithRawResponse,
    AsyncUniqueValuesResourceWithRawResponse,
    UniqueValuesResourceWithStreamingResponse,
    AsyncUniqueValuesResourceWithStreamingResponse,
)

__all__ = ["ColumnsResource", "AsyncColumnsResource"]


class ColumnsResource(SyncAPIResource):
    @cached_property
    def unique_values(self) -> UniqueValuesResource:
        return UniqueValuesResource(self._client)

    @cached_property
    def with_raw_response(self) -> ColumnsResourceWithRawResponse:
        return ColumnsResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> ColumnsResourceWithStreamingResponse:
        return ColumnsResourceWithStreamingResponse(self)


class AsyncColumnsResource(AsyncAPIResource):
    @cached_property
    def unique_values(self) -> AsyncUniqueValuesResource:
        return AsyncUniqueValuesResource(self._client)

    @cached_property
    def with_raw_response(self) -> AsyncColumnsResourceWithRawResponse:
        return AsyncColumnsResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncColumnsResourceWithStreamingResponse:
        return AsyncColumnsResourceWithStreamingResponse(self)


class ColumnsResourceWithRawResponse:
    def __init__(self, columns: ColumnsResource) -> None:
        self._columns = columns

    @cached_property
    def unique_values(self) -> UniqueValuesResourceWithRawResponse:
        return UniqueValuesResourceWithRawResponse(self._columns.unique_values)


class AsyncColumnsResourceWithRawResponse:
    def __init__(self, columns: AsyncColumnsResource) -> None:
        self._columns = columns

    @cached_property
    def unique_values(self) -> AsyncUniqueValuesResourceWithRawResponse:
        return AsyncUniqueValuesResourceWithRawResponse(self._columns.unique_values)


class ColumnsResourceWithStreamingResponse:
    def __init__(self, columns: ColumnsResource) -> None:
        self._columns = columns

    @cached_property
    def unique_values(self) -> UniqueValuesResourceWithStreamingResponse:
        return UniqueValuesResourceWithStreamingResponse(self._columns.unique_values)


class AsyncColumnsResourceWithStreamingResponse:
    def __init__(self, columns: AsyncColumnsResource) -> None:
        self._columns = columns

    @cached_property
    def unique_values(self) -> AsyncUniqueValuesResourceWithStreamingResponse:
        return AsyncUniqueValuesResourceWithStreamingResponse(
            self._columns.unique_values
        )
