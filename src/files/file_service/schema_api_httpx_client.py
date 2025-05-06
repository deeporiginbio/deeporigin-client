"""
Schema API HTTPX Client

A lightweight API client that integrates OpenAPI schema validation with httpx requests.
This client reads an OpenAPI schema file and uses it to validate requests and responses
against the schema definitions at runtime.

Key features:
- Direct JSON schema validation without generating Data Transfer Objects (DTOs)
- Reduces code bloat compared to traditional OpenAPI code generators
- Uses jsonschema library to validate API data against OpenAPI schema
- Handles content-type detection and appropriate parsing of responses

Usage:
    client = SchemaApiHttpxClient(
        base_url="https://api.example.com",
        schema_path="/path/to/openapi.yaml",
        token="optional-auth-token"
    )
    
    # Use with generated kwargs functions
    kwargs, path = api.get_some_endpoint_kwargs(param1="value")
    response = client.httpx_request_and_process((kwargs, path))
    
    # Check if response is valid according to schema
    if response.is_valid:
        data = response.parsed_json
    else:
        print(f"Validation error: {response.validation_error}")
"""

from dataclasses import dataclass
import json
from typing import Any, BinaryIO, Dict, Optional, Tuple

import httpx
import jsonschema
import yaml


@dataclass
class ApiResponse:
    """Simple wrapper for API responses that includes validation status"""

    status_code: int
    content: bytes
    headers: Dict[str, str]
    is_valid: bool = True
    parsed_json: Optional[Dict[str, Any]] = None
    validation_error: Optional[str] = None


class SchemaApiHttpxClient:
    """Base client for simplified API interactions"""

    _verify_ssl: str | bool
    _timeout: httpx.Timeout | float | None

    def __init__(
        self,
        base_url: str,
        schema_path: str,
        token: str = None,
        timeout: httpx.Timeout | float | None = 5.0,
        verify_ssl: bool | str = True,
    ):
        # Load the OpenAPI schema
        # schema_dir = os.path.dirname(os.path.abspath(__file__))
        with open(schema_path, "r") as f:
            self.schema = yaml.safe_load(f)

        # Setup HTTP client
        self.base_url = base_url
        self._verify_ssl = verify_ssl
        self._timeout = timeout
        self.client = httpx.Client(
            base_url=base_url, verify=self._verify_ssl, timeout=self._timeout
        )

        # Configure authentication if token provided
        if token:
            self.client.headers["Authorization"] = f"Bearer {token}"

    def get_request_schema(self, path: str, method: str) -> Dict:
        """Get request body schema for a specific endpoint"""
        try:
            path_obj = self.schema["paths"].get(path, {})
            method_obj = path_obj.get(method.lower(), {})
            return (
                method_obj.get("requestBody", {})
                .get("content", {})
                .get("application/json", {})
                .get("schema", {})
            )
        except (KeyError, TypeError):
            return {}

    def get_response_schema(self, path: str, method: str, status_code: str) -> Dict:
        """Get response schema for a specific endpoint and status code"""
        try:
            path_obj = self.schema["paths"].get(path, {})
            method_obj = path_obj.get(method.lower(), {})
            return (
                method_obj.get("responses", {})
                .get(status_code, {})
                .get("content", {})
                .get("application/json", {})
                .get("schema", {})
            )
        except (KeyError, TypeError):
            return {}

    def validate_request(
        self, data: Dict, path: str, method: str
    ) -> Tuple[bool, Optional[str]]:
        """Validate request data against schema"""
        schema = self.get_request_schema(path, method)
        if not schema:
            return True, None  # No schema to validate against

        try:
            jsonschema.validate(instance=data, schema=schema)
            return True, None
        except jsonschema.exceptions.ValidationError as e:
            return False, str(e)

    def validate_response(
        self, data: Dict, path: str, method: str, status_code: str
    ) -> Tuple[bool, Optional[str]]:
        """Validate response data against schema"""
        schema = self.get_response_schema(path, method, status_code)
        if not schema:
            return True, None  # No schema to validate against

        try:
            jsonschema.validate(instance=data, schema=schema)
            return True, None
        except jsonschema.exceptions.ValidationError as e:
            return False, str(e)

    def process_response(
        self, response: httpx.Response, path: str, method: str
    ) -> ApiResponse:
        """Process and validate response"""
        api_response = ApiResponse(
            status_code=response.status_code,
            content=response.content,
            headers=dict(response.headers),
        )

        # Try to parse JSON for responses with JSON content
        if response.headers.get("content-type", "").startswith("application/json"):
            try:
                api_response.parsed_json = response.json()

                # Validate against schema if status code is success
                if 200 <= response.status_code < 300:
                    is_valid, error = self.validate_response(
                        api_response.parsed_json,
                        path,
                        method,
                        str(response.status_code),
                    )
                    api_response.is_valid = is_valid
                    api_response.validation_error = error
            except json.JSONDecodeError:
                api_response.is_valid = False
                api_response.validation_error = "Invalid JSON response"

        return api_response

    def httpx_request_and_process(
        self, args_and_urlkey: tuple[dict, str]
    ) -> ApiResponse:
        """
        Perform an httpx request and process the response in one call
        """
        kwargs, url_key = args_and_urlkey
        response = self.client.request(**kwargs)
        return self.process_response(response, url_key, kwargs["method"])


# Utility classes for file handling


class FilePayload:
    """Wrapper for file payloads to be used in multipart requests"""

    def __init__(
        self,
        payload: BinaryIO,
        file_name: str = None,
        mime_type: str = "application/octet-stream",
    ):
        self.payload = payload
        self.file_name = file_name or "file"
        self.mime_type = mime_type

    def to_tuple(self) -> Tuple[str, Tuple[str, bytes, str]]:
        """Convert to tuple format expected by httpx for multipart/form-data"""
        # Httpx accepts BinaryIO for payload while streaming
        return ("file", (self.file_name, self.payload, self.mime_type))
