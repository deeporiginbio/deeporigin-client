"""module to test code in utils.types"""

import pytest
from beartype.roar import BeartypeCallHintParamViolation
from deeporigin.data_hub import filters
from pydantic import ValidationError

valid_operators = [
    "equals",
    "notEqual",
    "lessThan",
    "lessThanOrEqual",
    "greaterThan",
    "greaterThanOrEqual",
]


valid_numbers = [0, 1.2]

invalid_arguments = [
    ("foo", 100, "invalid"),
    (None, 100, "=="),  # invalid column_id (None)
    ("foo", "not_a_number", ">"),  # invalid value (not numeric)
    ("foo", 100, None),  # invalid operator (None)
    ("foo", 100, ""),  # invalid operator (empty string)
    (123, 100, "=="),  # invalid column_id (not a string)
]


@pytest.mark.parametrize("operator", valid_operators)
def test_numeric_condition_valid_operators(
    operator,
):
    filters.numeric_condition(
        column_id="foo",
        value=100,
        operator=operator,
    )


@pytest.mark.parametrize("number", valid_numbers)
def test_numeric_condition_valid_numbers(
    number,
):
    filters.numeric_condition(
        column_id="foo",
        value=number,
        operator=valid_operators[0],
    )


@pytest.mark.parametrize(
    "column_id, value, operator",
    invalid_arguments,
)
def test_numeric_condition_invalid(column_id, value, operator):
    with pytest.raises((ValidationError, BeartypeCallHintParamViolation)):
        # constructing this should raise an error
        filters.numeric_condition(
            column_id=column_id,
            value=value,
            operator=operator,
        )
