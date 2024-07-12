# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing import List

import httpx

from ..types import archive_file_archive_params
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

__all__ = ["ArchiveFilesResource", "AsyncArchiveFilesResource"]


class ArchiveFilesResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> ArchiveFilesResourceWithRawResponse:
        return ArchiveFilesResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> ArchiveFilesResourceWithStreamingResponse:
        return ArchiveFilesResourceWithStreamingResponse(self)

    def archive(
        self,
        *,
        file_ids: List[str],
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> object:
        """
        Archive files by their ids.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return self._post(
            "/ArchiveFiles",
            body=maybe_transform({"file_ids": file_ids}, archive_file_archive_params.ArchiveFileArchiveParams),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=object,
        )


class AsyncArchiveFilesResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncArchiveFilesResourceWithRawResponse:
        return AsyncArchiveFilesResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncArchiveFilesResourceWithStreamingResponse:
        return AsyncArchiveFilesResourceWithStreamingResponse(self)

    async def archive(
        self,
        *,
        file_ids: List[str],
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> object:
        """
        Archive files by their ids.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return await self._post(
            "/ArchiveFiles",
            body=await async_maybe_transform(
                {"file_ids": file_ids}, archive_file_archive_params.ArchiveFileArchiveParams
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=object,
        )


class ArchiveFilesResourceWithRawResponse:
    def __init__(self, archive_files: ArchiveFilesResource) -> None:
        self._archive_files = archive_files

        self.archive = to_raw_response_wrapper(
            archive_files.archive,
        )


class AsyncArchiveFilesResourceWithRawResponse:
    def __init__(self, archive_files: AsyncArchiveFilesResource) -> None:
        self._archive_files = archive_files

        self.archive = async_to_raw_response_wrapper(
            archive_files.archive,
        )


class ArchiveFilesResourceWithStreamingResponse:
    def __init__(self, archive_files: ArchiveFilesResource) -> None:
        self._archive_files = archive_files

        self.archive = to_streamed_response_wrapper(
            archive_files.archive,
        )


class AsyncArchiveFilesResourceWithStreamingResponse:
    def __init__(self, archive_files: AsyncArchiveFilesResource) -> None:
        self._archive_files = archive_files

        self.archive = async_to_streamed_response_wrapper(
            archive_files.archive,
        )
