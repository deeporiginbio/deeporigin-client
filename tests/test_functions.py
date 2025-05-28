"""this module contains tests for functions. These are meant to be run against a live instance"""

import pytest

from tests.utils import config  # noqa: F401


def test_molprops(config):  # noqa: F811
    if config["mock"]:
        pytest.skip("test skipped with mock client")

    from deeporigin.drug_discovery import Ligand

    ligand = Ligand.from_identifier("serotonin")

    props = ligand.admet_properties()

    assert isinstance(props, dict), "Expected a dictionary"
    assert "logP" in props, "Expected logP to be in the properties"
    assert "logD" in props, "Expected logD to be in the properties"
    assert "logS" in props, "Expected logS to be in the properties"
