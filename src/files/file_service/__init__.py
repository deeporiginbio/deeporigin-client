"""A client library for accessing File Service"""

from .schema_api_httpx_client import FilePayload, SchemaApiHttpxClient

__all__ = (
    "SchemaApiHttpxClient",
    "FilePayload",
)
