# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

import httpx

from ..types import import_row_create_params
from .._types import NOT_GIVEN, Body, Query, Headers, NotGiven
from .._utils import (
    maybe_transform,
    async_maybe_transform,
)
from .._compat import cached_property
from .._resource import SyncAPIResource, AsyncAPIResource
from .._response import (
    to_raw_response_wrapper,
    to_streamed_response_wrapper,
    async_to_raw_response_wrapper,
    async_to_streamed_response_wrapper,
)
from .._base_client import (
    make_request_options,
)
from ..types.import_row_create_response import ImportRowCreateResponse

__all__ = ["ImportRowsResource", "AsyncImportRowsResource"]


class ImportRowsResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> ImportRowsResourceWithRawResponse:
        return ImportRowsResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> ImportRowsResourceWithStreamingResponse:
        return ImportRowsResourceWithStreamingResponse(self)

    def create(
        self,
        *,
        database_id: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> ImportRowCreateResponse:
        """
        Creates new rows from CSV data.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return self._post(
            "/ImportRows",
            body=maybe_transform(
                {"database_id": database_id},
                import_row_create_params.ImportRowCreateParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers,
                extra_query=extra_query,
                extra_body=extra_body,
                timeout=timeout,
            ),
            cast_to=ImportRowCreateResponse,
        )


class AsyncImportRowsResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncImportRowsResourceWithRawResponse:
        return AsyncImportRowsResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncImportRowsResourceWithStreamingResponse:
        return AsyncImportRowsResourceWithStreamingResponse(self)

    async def create(
        self,
        *,
        database_id: str,
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> ImportRowCreateResponse:
        """
        Creates new rows from CSV data.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return await self._post(
            "/ImportRows",
            body=await async_maybe_transform(
                {"database_id": database_id},
                import_row_create_params.ImportRowCreateParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers,
                extra_query=extra_query,
                extra_body=extra_body,
                timeout=timeout,
            ),
            cast_to=ImportRowCreateResponse,
        )


class ImportRowsResourceWithRawResponse:
    def __init__(self, import_rows: ImportRowsResource) -> None:
        self._import_rows = import_rows

        self.create = to_raw_response_wrapper(
            import_rows.create,
        )


class AsyncImportRowsResourceWithRawResponse:
    def __init__(self, import_rows: AsyncImportRowsResource) -> None:
        self._import_rows = import_rows

        self.create = async_to_raw_response_wrapper(
            import_rows.create,
        )


class ImportRowsResourceWithStreamingResponse:
    def __init__(self, import_rows: ImportRowsResource) -> None:
        self._import_rows = import_rows

        self.create = to_streamed_response_wrapper(
            import_rows.create,
        )


class AsyncImportRowsResourceWithStreamingResponse:
    def __init__(self, import_rows: AsyncImportRowsResource) -> None:
        self._import_rows = import_rows

        self.create = async_to_streamed_response_wrapper(
            import_rows.create,
        )
