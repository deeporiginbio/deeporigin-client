# File generated from our OpenAPI spec by Stainless. See CONTRIBUTING.md for details.

from __future__ import annotations

from typing import Iterable

import httpx

from ..types import configure_column_select_option_configure_params
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
from ..types.configure_column_select_option_configure_response import ConfigureColumnSelectOptionConfigureResponse

__all__ = ["ConfigureColumnSelectOptionsResource", "AsyncConfigureColumnSelectOptionsResource"]


class ConfigureColumnSelectOptionsResource(SyncAPIResource):
    @cached_property
    def with_raw_response(self) -> ConfigureColumnSelectOptionsResourceWithRawResponse:
        return ConfigureColumnSelectOptionsResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> ConfigureColumnSelectOptionsResourceWithStreamingResponse:
        return ConfigureColumnSelectOptionsResourceWithStreamingResponse(self)

    def configure(
        self,
        *,
        column_id: str,
        option_configuration: Iterable[configure_column_select_option_configure_params.OptionConfiguration],
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> ConfigureColumnSelectOptionConfigureResponse:
        """Configure column select options.

        Supports both adding and removing options.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return self._post(
            "/ConfigureColumnSelectOptions",
            body=maybe_transform(
                {
                    "column_id": column_id,
                    "option_configuration": option_configuration,
                },
                configure_column_select_option_configure_params.ConfigureColumnSelectOptionConfigureParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=ConfigureColumnSelectOptionConfigureResponse,
        )


class AsyncConfigureColumnSelectOptionsResource(AsyncAPIResource):
    @cached_property
    def with_raw_response(self) -> AsyncConfigureColumnSelectOptionsResourceWithRawResponse:
        return AsyncConfigureColumnSelectOptionsResourceWithRawResponse(self)

    @cached_property
    def with_streaming_response(self) -> AsyncConfigureColumnSelectOptionsResourceWithStreamingResponse:
        return AsyncConfigureColumnSelectOptionsResourceWithStreamingResponse(self)

    async def configure(
        self,
        *,
        column_id: str,
        option_configuration: Iterable[configure_column_select_option_configure_params.OptionConfiguration],
        # Use the following arguments if you need to pass additional parameters to the API that aren't available via kwargs.
        # The extra values given here take precedence over values defined on the client or passed to this method.
        extra_headers: Headers | None = None,
        extra_query: Query | None = None,
        extra_body: Body | None = None,
        timeout: float | httpx.Timeout | None | NotGiven = NOT_GIVEN,
    ) -> ConfigureColumnSelectOptionConfigureResponse:
        """Configure column select options.

        Supports both adding and removing options.

        Args:
          extra_headers: Send extra headers

          extra_query: Add additional query parameters to the request

          extra_body: Add additional JSON properties to the request

          timeout: Override the client-level default timeout for this request, in seconds
        """
        return await self._post(
            "/ConfigureColumnSelectOptions",
            body=await async_maybe_transform(
                {
                    "column_id": column_id,
                    "option_configuration": option_configuration,
                },
                configure_column_select_option_configure_params.ConfigureColumnSelectOptionConfigureParams,
            ),
            options=make_request_options(
                extra_headers=extra_headers, extra_query=extra_query, extra_body=extra_body, timeout=timeout
            ),
            cast_to=ConfigureColumnSelectOptionConfigureResponse,
        )


class ConfigureColumnSelectOptionsResourceWithRawResponse:
    def __init__(self, configure_column_select_options: ConfigureColumnSelectOptionsResource) -> None:
        self._configure_column_select_options = configure_column_select_options

        self.configure = to_raw_response_wrapper(
            configure_column_select_options.configure,
        )


class AsyncConfigureColumnSelectOptionsResourceWithRawResponse:
    def __init__(self, configure_column_select_options: AsyncConfigureColumnSelectOptionsResource) -> None:
        self._configure_column_select_options = configure_column_select_options

        self.configure = async_to_raw_response_wrapper(
            configure_column_select_options.configure,
        )


class ConfigureColumnSelectOptionsResourceWithStreamingResponse:
    def __init__(self, configure_column_select_options: ConfigureColumnSelectOptionsResource) -> None:
        self._configure_column_select_options = configure_column_select_options

        self.configure = to_streamed_response_wrapper(
            configure_column_select_options.configure,
        )


class AsyncConfigureColumnSelectOptionsResourceWithStreamingResponse:
    def __init__(self, configure_column_select_options: AsyncConfigureColumnSelectOptionsResource) -> None:
        self._configure_column_select_options = configure_column_select_options

        self.configure = async_to_streamed_response_wrapper(
            configure_column_select_options.configure,
        )
