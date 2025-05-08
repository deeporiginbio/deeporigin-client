# Schema API Generator
#
# This script generates functions based on an OpenAPI schema (.yaml), to support
# the use of lightweight SchemaApiHttpxClient. It produces simple Python functions
# that correspond to each API endpoint defined in the schema.
# We avoid the bloat of Data Transfer Objects, created by other OpenAPI generators,
# expecting simple json validation to do the job instead.
#
# Each generated function:
# - Takes parameters matching the API endpoint parameters
# - Constructs a proper request configuration dictionary for httpx
# - Returns both the request parameters and original path template (for schema lookup)
#
# The generated functions follow the naming pattern _get_operation_id_kwargs() and
# return tuple[dict[str, any], str] where:
# - First element: Dictionary with request kwargs for httpx
# - Second element: Original path template string for schema validation

import argparse
import os
import re

import yaml


def camel_to_snake(name: str) -> str:
    """Convert camelCase to snake_case"""
    name = re.sub("(.)([A-Z][a-z]+)", r"\1_\2", name)
    return re.sub("([a-z0-9])([A-Z])", r"\1_\2", name).lower()


# Track generated function names to handle duplicates
_generated_function_names = {}


def get_unique_function_name(operation_id: str, path: str) -> str:
    """Generate a unique function name for an operation"""
    # Convert operation_id to snake case
    base_name = camel_to_snake(operation_id)
    full_name = f"_get_{base_name}_kwargs"
    
    # Get the path prefix and convert to valid Python identifier
    path_prefix = path.split("/")[1]  # e.g., "file-api/shared-data" -> "file-api"
    path_prefix = path_prefix.replace("-", "_")  # Convert hyphens to underscores
    
    # Check if this function name has been used before
    if full_name in _generated_function_names:
        # Append the path prefix to make it unique
        unique_name = f"_get_{base_name}_{path_prefix}_kwargs"
        _generated_function_names[unique_name] = True
        return unique_name
    else:
        _generated_function_names[full_name] = True
        return full_name


def format_path_as_name(path: str, method: str) -> str:
    """Convert path and method to function name"""
    # Replace path params like {org_id} with "by_org_id"
    path = re.sub(r"\{([^}]+)\}", lambda m: f"by_{m.group(1)}", path)

    # Replace slashes and hyphens with underscores
    path = path.replace("/", "_").replace("-", "_")

    # Remove leading and trailing underscores
    return path.strip("_")


def sanitize_param_name(name: str) -> str:
    """Convert parameter name to valid Python identifier by replacing dashes with underscores."""
    return name.replace("-", "_")


def generate_api_file(schema_path: str, output_path: str) -> None:
    """Generate API client file from OpenAPI schema"""
    # Reset the function name tracking
    global _generated_function_names
    _generated_function_names = {}

    # Load OpenAPI schema
    with open(schema_path, "r") as f:
        schema = yaml.safe_load(f)

    # Start building output file content
    output = [
        "# Generated API client - DO NOT EDIT MANUALLY",
        "# This file was auto-generated from an OpenAPI schema",
        "",
        "import httpx",
        "from .schema_api_httpx_client import FilePayload",
        "",
        "",
    ]

    # Process each path and method in the schema
    for path, path_item in schema.get("paths", {}).items():
        for method, operation in path_item.items():
            if method not in ["get", "post", "put", "delete", "head"]:
                continue

            operation_id = operation.get(
                "operationId", format_path_as_name(path, method)
            )
            # Get unique function name
            function_name = get_unique_function_name(operation_id, path)

            # Get parameters from the operation
            parameters = operation.get("parameters", [])

            # Get request body schema if present
            request_body = operation.get("requestBody", {})
            has_multipart = False
            if request_body:
                content_types = request_body.get("content", {})
                has_multipart = "multipart/form-data" in content_types

            # Build function
            func_lines = [f"def {function_name}("]

            # Add path parameters
            path_params = []
            for param in parameters:
                if param.get("in") == "path":
                    param_name = param.get("name")
                    path_params.append(param_name)
                    func_lines.append(f"    {param_name}: str,")

            # Add query parameters
            query_params = []
            for param in parameters:
                if param.get("in") == "query":
                    param_name = param.get("name")
                    python_param_name = sanitize_param_name(param_name)
                    query_params.append((param_name, python_param_name))
                    param_type = "str"  # Default type
                    if param.get("schema", {}).get("type") == "number":
                        param_type = "float"
                    func_lines.append(
                        f"    {python_param_name}: {param_type} | None = None,"
                    )

            # Add request body if present
            if has_multipart:
                func_lines.append("    file: FilePayload | None = None,")
            elif request_body:
                func_lines.append("    body: dict[str, any] | None = None,")

            # Close function signature with tuple return type
            func_lines.append(") -> tuple[dict[str, any], str]:")

            # Add short function docstring
            func_lines.append(f'    """Get request kwargs for {operation_id}"""')

            # Add function body
            func_lines.append("    _kwargs: dict[str, any] = {")
            func_lines.append(f'        "method": "{method}",')

            # Format the URL with path parameters
            #url = path
            if path_params:
                func_lines.append(
                    f'        "url": f"{path}".format({", ".join(f"{p}={p}" for p in path_params)}),'
                )
            else:
                func_lines.append(f'        "url": "{path}",')

            func_lines.append("    }")

            # Add query parameters
            if query_params:
                func_lines.append("    params: dict[str, any] = {}")
                for orig_name, python_name in query_params:
                    func_lines.append(f"    if {python_name} is not None:")
                    func_lines.append(f'        params["{orig_name}"] = {python_name}')
                func_lines.append("    if params:")
                func_lines.append('        _kwargs["params"] = params')

            # Add request body
            if has_multipart:
                func_lines.append("    if file is not None:")
                func_lines.append(
                    '        _kwargs["files"] = {"file": file.to_tuple()[1]}'
                )
            elif request_body:
                func_lines.append("    if body is not None:")
                func_lines.append('        _kwargs["json"] = body')
                func_lines.append(
                    '        _kwargs["headers"] = {"Content-Type": "application/json"}'
                )

            # Return kwargs and original path as a tuple
            func_lines.append(f'    return _kwargs, "{path}"')
            func_lines.append("")

            # Add function to output
            output.extend(func_lines)

    # Write output file
    with open(output_path, "w") as f:
        f.write("\n".join(output))

    print(f"Generated API client at {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Generate API client from OpenAPI schema"
    )
    parser.add_argument(
        "--schema",
        "-s",
        help="Path to OpenAPI schema YAML file (default: gen_file_openapi.yaml in current directory)",
    )
    args = parser.parse_args()

    script_dir = os.path.dirname(os.path.abspath(__file__))

    # Use command line argument if provided, otherwise look for default file in current directory
    if args.schema:
        schema_path = args.schema
    else:
        # First check if default file exists in current directory
        default_file = "gen_file_openapi.yaml"
        if os.path.exists(default_file):
            schema_path = default_file
        else:
            # Fall back to the original relative path
            schema_path = os.path.join(
                script_dir, "../../scripts/gen_file_openapi.yaml"
            )

    output_path = os.path.join(script_dir, "generated_api.py")

    print(f"Using schema file: {schema_path}")
    generate_api_file(schema_path, output_path)
