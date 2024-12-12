"""module to create filters and conditions to filter rows in databases on the data hub"""

from typing import Any, Literal

from beartype import beartype
from box import Box

FilterType = Literal["number", "text", "boolean"]
Operator = Literal[
    "equals",
    "notEqual",
    "lessThan",
    "lessThanOrEqual",
    "greaterThan",
    "greaterThanOrEqual",
    "contains",
    "notContains",
    "startsWith",
    "endsWith",
    "substructure",
    "in",
    "notIn",
    "isNull",
    "isNotNull",
]


@beartype
def filter(
    *,
    filter_type: FilterType,
    operator: Operator,
    filter_value: Any,
    column_id: str,
):
    """function to create filters and conditions to filter rows in databases on the data hub

    Args:
        filter_type (FilterType): type of filter. can be one of number, text, boolean
        operator (Operator): operator to apply
        filter_value (Any): value of filter
        column_id (str): id or name of column


    """
    data = dict(
        column_id=column_id,
        filter_type=filter_type,
        filter_value=filter_value,
        operator=operator,
    )

    return Box(data)
