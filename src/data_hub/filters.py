"""module to create filters and conditions to filter rows in databases on the data hub"""

from beartype import beartype
from deeporigin_data.types.client_list_database_rows_params import (
    FilterRowFilterNumber,
    FilterRowFilterText,
)
from pydantic import TypeAdapter


@beartype
def numeric_condition(
    *,
    column_id: str,
    value: float | int,
    operator: str,
) -> dict:
    """function to create a valid numeric filter


    Args:
        column_id (str): id of the column
        value (float | int): value to filter
        operator (str): operator to use.

    Returns:
        A dictionary containing a valid condition that can be used to filter database rows
    """

    data = dict(
        column_id=column_id,
        filter_type="number",
        filter_value=value,
        operator=operator,
    )

    # use type checks from SDK definition, which
    # ultimately derives from the openAPI spec
    adapter = TypeAdapter(FilterRowFilterNumber)
    adapter.validate_python(data)

    return data


@beartype
def text_condition(
    *,
    column_id: str,
    value: str,
    operator: str,
) -> dict:
    """function to create a valid text filter"""
    data = dict(
        column_id=column_id,
        filter_type="text",
        filter_value=value,
        operator=operator,
    )

    # use type checks from SDK definition, which
    # ultimately derives from the openAPI spec
    adapter = TypeAdapter(FilterRowFilterText)
    adapter.validate_python(data)

    return data
