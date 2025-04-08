"""Visualization utilities for drug discovery using Jupyter notebooks.

This module provides decorators and utilities for visualizing molecular structures
and other drug discovery related data in Jupyter notebooks using the DeepOrigin Molstar viewer.
"""

from deeporigin_molstar import JupyterViewer


def jupyter_visualization(func):
    """Decorator for converting HTML visualization output to Jupyter notebook display.

    This decorator wraps functions that generate HTML visualizations and converts
    their output to be properly displayed in Jupyter notebooks using the DeepOrigin
    Molstar viewer.

    Args:
        func (callable): A function that returns HTML visualization content.

    Returns:
        callable: A wrapped function that returns a JupyterViewer visualization.

    Example:
        @jupyter_visualization
        def generate_molecule_view(molecule):
            # Generate HTML visualization
            return html_content
    """

    def wrapper(*args, **kwargs):
        html_visualization = func(*args, **kwargs)
        return JupyterViewer.visualize(html_visualization)

    return wrapper
