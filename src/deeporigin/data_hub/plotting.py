"""module to make plots from a Deep Origin dataframe"""

from beartype import beartype
from beartype.typing import Optional
from bokeh.io import show
from bokeh.layouts import Spacer, column, row
from bokeh.models import (
    Button,
    ColumnDataSource,
    CustomJS,
    DataTable,
    HoverTool,
    LassoSelectTool,
    Select,
    TableColumn,
)
from bokeh.palettes import Category10
from bokeh.plotting import figure
import pandas as pd

from deeporigin.data_hub.dataframe import DataFrame
from deeporigin.exceptions import DeepOriginException


@beartype
def scatter(
    df: DataFrame,
    x: Optional[str] = None,
    y: Optional[str] = None,
    size: Optional[str] = None,
    label_column: Optional[str] = None,
    js_code: Optional[dict] = None,
):
    """function to make a scatter plot from a Deep Origin dataframe, with support for interactivity

    Args:
        df (DataFrame): Deep Origin dataframe
        x (Optional[str], optional): name of column to use for x axis. Defaults to None.
        y (Optional[str], optional): name of column to use for y axis. Defaults to None.
        size (Optional[str], optional): name of column to use for size. Defaults to None.
        hover_callback_code: Optional[str], optional): JavaScript callback to use for hover tool. Defaults to None.

    `df` should be a Deep Origin dataframe with at least two numeric columns.

    """

    # constants for this plot
    figure_width = 500
    select_width = int(figure_width * 0.3)
    if js_code is None:
        js_code = _read_js_code()
    default_label = "None"
    default_color = "blue"

    cols = df.attrs["metadata"]["cols"]
    numeric_cols = [col["name"] for col in cols if col["type"] in ["float", "integer"]]

    # make placeholder lists for legends and colors
    legend_labels = [default_label for _ in range(len(df))]
    colors = [default_color for _ in range(len(df))]

    if label_column is not None:
        select_cols = [
            col for col in df.attrs["metadata"]["cols"] if col["type"] == "select"
        ]
        if len(select_cols) == 0:
            raise DeepOriginException("No select columns detected in the dataframe.")
        select_cols = [col["name"] for col in select_cols]
        if label_column not in select_cols:
            raise DeepOriginException(
                f"Column '{label_column}' is not a select column in the dataframe."
            )

        labels = df[label_column].cat.categories.tolist()

        if len(labels) == 0:
            raise DeepOriginException(
                f"The column {label_column} has no options configured, and therefore cannot be used to label data."
            )

        if len(labels) > 10:
            raise DeepOriginException(
                f"The column {label_column} has too many ({len(labels)}) options configured, and therefore cannot be used to label data."
            )

        # make color map
        cat10_colors = Category10
        cat10_colors[1] = cat10_colors[3][:1]
        cat10_colors[2] = cat10_colors[3][:2]
        cat10_colors = cat10_colors[len(labels)]
        color_map = {
            label: color for label, color in zip(labels, cat10_colors, strict=False)
        }
        color_map[default_label] = default_color

        legend_labels = [
            default_label if pd.isna(x) else x for x in df[label_column].tolist()
        ]

        colors = [color_map[label] for label in legend_labels]

    if len(numeric_cols) < 2:
        raise DeepOriginException(
            "DataFrame must contain at least two numeric columns."
        )

    if x is None:
        x = numeric_cols[0]
    if y is None:
        y = numeric_cols[1]
    if size is None:
        size = numeric_cols[0]

    # normalize sizes. this should match what's in
    # axes_callback.js
    sizes = list(df[size])
    min_size = min(sizes)
    max_size = max(sizes)
    sizes = [2 + 15 * (value - min_size) / (max_size - min_size) for value in sizes]

    # CDS for scatter data
    data = dict(
        x=list(df[x]),
        y=list(df[y]),
        size=sizes,
        id=df.index,
        colors=colors,
        legend_labels=legend_labels,
    )

    # also add every column to the data
    for col in df.columns:
        data[col] = list(df[col])

    # CDS for scatter plot
    scatter_source = ColumnDataSource(data)

    # we want to make a table with two columns that will
    # show the currently hovered point
    table_source = ColumnDataSource(data={"Column Name": [], "Value": []})
    columns = [
        TableColumn(field="Column Name", title="Column Name"),
        TableColumn(field="Value", title="Value"),
    ]
    data_table = DataTable(
        source=table_source,
        columns=columns,
        width=250,
        height=figure_width,
        index_position=None,
    )

    # CDS for marker
    marker_source = ColumnDataSource(_first_element_in_dict(data))

    # CDS to store lasso selection data
    lasso_selection_source = ColumnDataSource(dict(ids=[]))

    # Create the scatter plot figure
    p = figure(
        width=figure_width,
        height=figure_width,
        title="Deep Origin / " + df.attrs["metadata"]["name"],
    )
    p.toolbar.logo = None

    # make scatter plot of data
    scatter_glyphs = p.scatter(
        x="x",
        y="y",
        source=scatter_source,
        size="size",
        color="colors",
        alpha=0.6,
        legend_field="legend_labels",
    )

    # draw the marker
    p.scatter(
        x="x",
        y="y",
        size=10,
        fill_color=None,
        color="red",
        source=marker_source,
        line_width=3,
    )

    # Set initial axis labels
    p.xaxis.axis_label = x
    p.yaxis.axis_label = y

    # Create dropdown selectors for X and Y axes
    x_select = Select(
        title="X-Axis",
        value=x,
        options=numeric_cols,
        width=select_width,
    )
    y_select = Select(
        title="Y-Axis",
        value=y,
        options=numeric_cols,
        width=select_width,
    )
    size_select = Select(
        title="Size",
        value=size,
        options=numeric_cols,
        width=select_width,
    )

    # create dropdown selectors for label column
    if label_column:
        label_select = Select(
            title="Label",
            value=labels[0],
            options=labels,
            width=select_width,
        )

    # JavaScript callback to update data, axis labels, and point sizes on select change

    axes_callback = CustomJS(
        args=dict(
            scatter_source=scatter_source,
            x_select=x_select,
            y_select=y_select,
            size_select=size_select,
            df=df.reset_index().to_dict("list"),
            x_axis=p.xaxis[0],
            y_axis=p.yaxis[0],
        ),
        code=js_code["axes_callback"],
    )

    # JS code, will run in browser
    # this updates the value of the slider to the currently
    # hovered point

    hover_callback = CustomJS(
        code=js_code["hover_callback"],
        args=dict(
            marker_source=marker_source,
            scatter_source=scatter_source,
            table_source=table_source,
        ),
    )

    # custom hover tool to ignore marker glyphs
    # see this for an explanation of why we do this
    # https://discourse.bokeh.org/t/deactivate-hovertool-for-specific-glyphs/9931/2
    hvr = HoverTool(
        tooltips=None,
        callback=hover_callback,
    )
    hvr.renderers = [scatter_glyphs]
    p.add_tools(hvr)

    # Attach the callback to the select widgets
    x_select.js_on_change("value", axes_callback)
    y_select.js_on_change("value", axes_callback)
    size_select.js_on_change("value", axes_callback)

    lasso_callback = CustomJS(
        args=dict(
            scatter_source=scatter_source,
            lasso_selection_source=lasso_selection_source,
        ),
        code=js_code["lasso_callback"],
    )

    lasso_tool = LassoSelectTool()
    p.add_tools(lasso_tool)
    scatter_source.selected.js_on_change("indices", lasso_callback)

    # legend
    p.legend.title = "Label"

    # Button to access selected data from selected_source
    select_row = row(x_select, y_select, size_select)
    if label_column:
        label_button = Button(label="Add")
        button_callback = CustomJS(
            args=dict(
                lasso_selection_source=lasso_selection_source,
                label_select=label_select,
                scatter_source=scatter_source,
                color_map=color_map,
                label_column=label_column,
            ),
            code=js_code["button_callback"],
        )
        label_button.js_on_click(button_callback)
        select_row.children.append(Spacer(width=30))
        select_row.children.append(label_select)
        select_row.children.append(
            column(
                Spacer(width=10, height=20),
                label_button,
            )
        )

    layout = column(
        select_row,
        row(p, data_table),
    )
    show(layout)


@beartype
def _first_element_in_dict(data: dict) -> dict:
    """utility function that only includes the first element
    of lists in a dict, useful for building a marker-based
    CDS"""

    out_data = dict()
    for key in data.keys():
        out_data[key] = [data[key][0]]

    return out_data


@beartype
def _read_js_code() -> dict:
    """utility function to read JS code"""

    files = _list_files("js")
    js_code = dict()

    for key in files.keys():
        with open(files[key], "r") as f:
            js_code[key] = f.read()

    return js_code


@beartype
def _list_files(ext: str = "js") -> dict:
    """utility function that lists all files by extension in the package data"""

    import importlib.resources

    files = []
    with importlib.resources.files("deeporigin") as package_files:
        for file in package_files.rglob("*." + ext):
            files.append(file)

    return {path.stem: path for path in files}
