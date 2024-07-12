# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

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
from ...types.rows import back_reference_list_params
from ..._base_client import (
    make_request_options,
)
from ...types.rows.back_reference_list_response import BackReferenceListResponse

__all__ = ["BackReferencesResource", "AsyncBackReferencesResource"]


class BackReferencesResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> BackReferencesResourceWithRawResponse:
        return BackReferencesResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> BackReferencesResourceWithStreamingResponse:
        return BackReferencesResourceWithStreamingResponse(self)

    def list(
        self,
        *,
        row_id: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> BackReferenceListResponse:
        """
        Finds all the places a row is referenced.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return self._post(
            "/ListRowBackReferences",
            body=maybe_transform(
                {"row_id": row_id}, back_reference_list_params.BackReferenceListParams
            ),
            options=make_request_options(
                extra_headers=extra_headers,
                extra_query=extra_query,
                extra_body=extra_body,
                timeout=timeout,
            ),
            cast_to=BackReferenceListResponse,
        )


class AsyncBackReferencesResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncBackReferencesResourceWithRawResponse:
        return AsyncBackReferencesResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(
        self,
    ) -> AsyncBackReferencesResourceWithStreamingResponse:
        return AsyncBackReferencesResourceWithStreamingResponse(self)

    async def list(
        self,
        *,
        row_id: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> BackReferenceListResponse:
        """
        Finds all the places a row is referenced.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return await self._post(
            "/ListRowBackReferences",
            body=await async_maybe_transform(
                {"row_id": row_id}, back_reference_list_params.BackReferenceListParams
            ),
            options=make_request_options(
                extra_headers=extra_headers,
                extra_query=extra_query,
                extra_body=extra_body,
                timeout=timeout,
            ),
            cast_to=BackReferenceListResponse,
        )


class BackReferencesResourceWithRawResponse:
    def __init__(self, back_references: BackReferencesResource) -> None:
        self._back_references = back_references

        self.list = to_raw_response_wrapper(
            back_references.list,
        )


class AsyncBackReferencesResourceWithRawResponse:
    def __init__(self, back_references: AsyncBackReferencesResource) -> None:
        self._back_references = back_references

        self.list = async_to_raw_response_wrapper(
            back_references.list,
        )


class BackReferencesResourceWithStreamingResponse:
    def __init__(self, back_references: BackReferencesResource) -> None:
        self._back_references = back_references

        self.list = to_streamed_response_wrapper(
            back_references.list,
        )


class AsyncBackReferencesResourceWithStreamingResponse:
    def __init__(self, back_references: AsyncBackReferencesResource) -> None:
        self._back_references = back_references

        self.list = async_to_streamed_response_wrapper(
            back_references.list,
        )
