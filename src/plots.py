"""this module contains functions for plotting"""

import math
from typing import Optional, Sequence

from bokeh.io import show
from bokeh.models import (
    BasicTicker,
    ColorBar,
    ColumnDataSource,
    HoverTool,
    LinearColorMapper,
    PrintfTickFormatter,
)
from bokeh.palettes import Viridis256
from bokeh.plotting import figure
import numpy as np


def plot_heatmap(
    values: np.ndarray,
    *,
    labels: Optional[Sequence[str]] = None,
    title: str = "",
    palette=Viridis256,
    size: int = 700,
    show_values_on_hover: bool = True,
    clim: Optional[tuple[float, float]] = None,
):
    """
    Visualize a square matrix (NxN) as a Bokeh heatmap.

    Parameters
    ----------
    values : np.ndarray
        Square NxN matrix. NaNs are allowed.
    labels : list[str], optional
        Row/column labels. If None, uses "0..N-1".
    title : str
        Plot title.
    palette : sequence of colors
        Bokeh palette for the heatmap.
    size : int
        Figure size in pixels (width = height for square matrix).
    show_values_on_hover : bool
        If True, shows (row, col, value) tooltips when hovering cells.
    clim : tuple[float, float], optional
        Color limits as (vmin, vmax). If None, automatically computed from data.
        Useful for consistent color scaling across multiple plots.

    Returns
    -------
    bokeh.plotting.Figure
    """
    # --- Normalize inputs ---
    mat = np.asarray(values, dtype=float)
    if labels is None:
        labels = [str(i) for i in range(mat.shape[0])]
    else:
        labels = list(map(str, labels))

    if mat.ndim != 2 or mat.shape[0] != mat.shape[1]:
        raise ValueError("R must be a square NxN matrix.")

    n = mat.shape[0]
    if len(labels) != n:
        raise ValueError("Length of `labels` must match the matrix dimension.")

    # --- Compute color scale range, ignoring NaNs ---
    if clim is not None:
        vmin, vmax = clim
    else:
        finite_vals = mat[np.isfinite(mat)]
        if finite_vals.size == 0:
            vmin, vmax = 0.0, 1.0
        else:
            vmin, vmax = float(np.nanmin(finite_vals)), float(np.nanmax(finite_vals))
            if math.isclose(vmin, vmax):
                # Avoid degenerate color scales
                delta = 1e-6 if vmin == 0 else abs(vmin) * 1e-6
                vmin, vmax = vmin - delta, vmax + delta

    # --- Build rect grid data (categorical axes so labels show nicely) ---
    xs, ys, vals, ii, jj = [], [], [], [], []
    for i in range(n):
        for j in range(n):
            xs.append(labels[j])  # x is column
            ys.append(labels[n - 1 - i])  # y reversed so (0,0) is top-left visually
            vals.append(mat[i, j])
            ii.append(i)
            jj.append(j)

    source = ColumnDataSource(
        {
            "x": xs,
            "y": ys,
            "value": vals,
            "i": ii,
            "j": jj,
            "value_str": [("NA" if not np.isfinite(v) else f"{v:.4f}") for v in vals],
        }
    )

    # --- Color mapper ---
    mapper = LinearColorMapper(
        palette=palette,
        low=vmin,
        high=vmax,
        nan_color="#dddddd",
    )

    # --- Create and configure figure ---
    p = figure(
        title=title,
        x_range=labels,
        y_range=list(reversed(labels)),
        x_axis_location="above",
        tools="pan,wheel_zoom,box_zoom,reset,save",
        toolbar_location="right",
        width=size,
        height=size,
        tooltips=None,
        match_aspect=True,
    )

    # Add cells
    p.rect(
        x="x",
        y="y",
        width=1,
        height=1,
        source=source,
        line_color=None,
        fill_color={"field": "value", "transform": mapper},
    )

    if show_values_on_hover:
        hover = HoverTool(
            tooltips=[
                ("row (i)", "@i"),
                ("col (j)", "@j"),
                ("label row", "@y"),
                ("label col", "@x"),
                ("RMSD", "@value_str"),
            ]
        )
        p.add_tools(hover)

    color_bar = ColorBar(
        color_mapper=mapper,
        location=(0, 0),
        ticker=BasicTicker(desired_num_ticks=8),
        formatter=PrintfTickFormatter(format="%.3f"),
        label_standoff=8,
    )
    p.add_layout(color_bar, "right")

    p.axis.major_label_text_font_size = "9pt"
    p.xaxis.major_label_orientation = 0.9
    p.grid.visible = False

    show(p)
