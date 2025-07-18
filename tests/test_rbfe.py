"""this module contains tests for the rbfe module"""

import pytest

from deeporigin.drug_discovery import BRD_DATA_DIR, Complex
from deeporigin.exceptions import DeepOriginException


def test_network_should_exist():
    """test that the network exists"""

    sim = Complex.from_dir(BRD_DATA_DIR)

    with pytest.raises(
        DeepOriginException,
        match="Network not mapped yet.",
    ):
        sim.rbfe.run_network()
