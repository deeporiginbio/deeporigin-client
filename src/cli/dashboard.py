import subprocess
import sys


def run():
    """Run the Streamlit app in the current virtual environment."""
    subprocess.run([sys.executable, "-m", "streamlit", "run", "src/tools/dashboard.py"])
