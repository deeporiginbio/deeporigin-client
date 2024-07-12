# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing import Iterable

import httpx

from ..._types import NOT_GIVEN, Body, Query, Headers, NotGiven
from ..._utils import (
    maybe_transform,
    async_maybe_transform,
)
from ..._compat import cached_property
from ..._resource import SyncAPIResource, AsyncAPIResource
from ..._response import (
    to_raw_response_wrapper,
    to_streamed_response_wrapper,
    async_to_raw_response_wrapper,
    async_to_streamed_response_wrapper,
)
from ...types.rows import hierarchy_list_params
from ..._base_client import (
    make_request_options,
)
from ...types.rows.hierarchy_list_response import HierarchyListResponse

__all__ = ["HierarchyResource", "AsyncHierarchyResource"]


class HierarchyResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> HierarchyResourceWithRawResponse:
        return HierarchyResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> HierarchyResourceWithStreamingResponse:
        return HierarchyResourceWithStreamingResponse(self)

    def list(
        self,
        *,
        filters: Iterable[hierarchy_list_params.Filter],
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> HierarchyListResponse:
        """
        Lists rows at a given depth in the hierarchy.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return self._post(
            "/ListRows",
            body=maybe_transform({"filters": filters}, hierarchy_list_params.HierarchyListParams),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=HierarchyListResponse,
        )


class AsyncHierarchyResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncHierarchyResourceWithRawResponse:
        return AsyncHierarchyResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncHierarchyResourceWithStreamingResponse:
        return AsyncHierarchyResourceWithStreamingResponse(self)

    async def list(
        self,
        *,
        filters: Iterable[hierarchy_list_params.Filter],
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> HierarchyListResponse:
        """
        Lists rows at a given depth in the hierarchy.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return await self._post(
            "/ListRows",
            body=await async_maybe_transform({"filters": filters}, hierarchy_list_params.HierarchyListParams),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=HierarchyListResponse,
        )


class HierarchyResourceWithRawResponse:
    def __init__(self, hierarchy: HierarchyResource) -> None:
        self._hierarchy = hierarchy

        self.list = to_raw_response_wrapper(
            hierarchy.list,
        )


class AsyncHierarchyResourceWithRawResponse:
    def __init__(self, hierarchy: AsyncHierarchyResource) -> None:
        self._hierarchy = hierarchy

        self.list = async_to_raw_response_wrapper(
            hierarchy.list,
        )


class HierarchyResourceWithStreamingResponse:
    def __init__(self, hierarchy: HierarchyResource) -> None:
        self._hierarchy = hierarchy

        self.list = to_streamed_response_wrapper(
            hierarchy.list,
        )


class AsyncHierarchyResourceWithStreamingResponse:
    def __init__(self, hierarchy: AsyncHierarchyResource) -> None:
        self._hierarchy = hierarchy

        self.list = async_to_streamed_response_wrapper(
            hierarchy.list,
        )
