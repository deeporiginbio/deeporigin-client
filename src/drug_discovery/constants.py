"""Constants and configuration values for drug discovery related functionality.

This module contains various constants used throughout the drug discovery pipeline,
including lists of metal atoms and file paths for state management.
"""

from pathlib import Path

METALS = ["MG", "CA", "ZN", "MN", "CU", "FE", "NA", "K", "HG", "CO", "U", "CD", "NI"]

# File paths
STATE_DUMP_PATH = Path.home() / ".deeporigin" / "state_dump.pdb"
