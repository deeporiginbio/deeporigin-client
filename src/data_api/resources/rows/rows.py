# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from ..._compat import cached_property
from .hierarchy import (
    HierarchyResource,
    AsyncHierarchyResource,
    HierarchyResourceWithRawResponse,
    AsyncHierarchyResourceWithRawResponse,
    HierarchyResourceWithStreamingResponse,
    AsyncHierarchyResourceWithStreamingResponse,
)
from ..._resource import SyncAPIResource, AsyncAPIResource
from .back_references import (
    BackReferencesResource,
    AsyncBackReferencesResource,
    BackReferencesResourceWithRawResponse,
    AsyncBackReferencesResourceWithRawResponse,
    BackReferencesResourceWithStreamingResponse,
    AsyncBackReferencesResourceWithStreamingResponse,
)

__all__ = ["RowsResource", "AsyncRowsResource"]


class RowsResource(SyncAPIResource):
    @cached_property
    def back_references(self) -> BackReferencesResource:
        return BackReferencesResource(self._client)

    @cached_property
    def hierarchy(self) -> HierarchyResource:
        return HierarchyResource(self._client)

    @cached_property
    def with_raw_response(self) -> RowsResourceWithRawResponse:
        return RowsResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> RowsResourceWithStreamingResponse:
        return RowsResourceWithStreamingResponse(self)


class AsyncRowsResource(AsyncAPIResource):
    @cached_property
    def back_references(self) -> AsyncBackReferencesResource:
        return AsyncBackReferencesResource(self._client)

    @cached_property
    def hierarchy(self) -> AsyncHierarchyResource:
        return AsyncHierarchyResource(self._client)

    @cached_property
    def with_raw_response(self) -> AsyncRowsResourceWithRawResponse:
        return AsyncRowsResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncRowsResourceWithStreamingResponse:
        return AsyncRowsResourceWithStreamingResponse(self)


class RowsResourceWithRawResponse:
    def __init__(self, rows: RowsResource) -> None:
        self._rows = rows

    @cached_property
    def back_references(self) -> BackReferencesResourceWithRawResponse:
        return BackReferencesResourceWithRawResponse(self._rows.back_references)

    @cached_property
    def hierarchy(self) -> HierarchyResourceWithRawResponse:
        return HierarchyResourceWithRawResponse(self._rows.hierarchy)


class AsyncRowsResourceWithRawResponse:
    def __init__(self, rows: AsyncRowsResource) -> None:
        self._rows = rows

    @cached_property
    def back_references(self) -> AsyncBackReferencesResourceWithRawResponse:
        return AsyncBackReferencesResourceWithRawResponse(self._rows.back_references)

    @cached_property
    def hierarchy(self) -> AsyncHierarchyResourceWithRawResponse:
        return AsyncHierarchyResourceWithRawResponse(self._rows.hierarchy)


class RowsResourceWithStreamingResponse:
    def __init__(self, rows: RowsResource) -> None:
        self._rows = rows

    @cached_property
    def back_references(self) -> BackReferencesResourceWithStreamingResponse:
        return BackReferencesResourceWithStreamingResponse(self._rows.back_references)

    @cached_property
    def hierarchy(self) -> HierarchyResourceWithStreamingResponse:
        return HierarchyResourceWithStreamingResponse(self._rows.hierarchy)


class AsyncRowsResourceWithStreamingResponse:
    def __init__(self, rows: AsyncRowsResource) -> None:
        self._rows = rows

    @cached_property
    def back_references(self) -> AsyncBackReferencesResourceWithStreamingResponse:
        return AsyncBackReferencesResourceWithStreamingResponse(self._rows.back_references)

    @cached_property
    def hierarchy(self) -> AsyncHierarchyResourceWithStreamingResponse:
        return AsyncHierarchyResourceWithStreamingResponse(self._rows.hierarchy)
