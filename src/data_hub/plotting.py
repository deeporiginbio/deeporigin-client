"""module to make plots from a Deep Origin dataframe"""

from beartype import beartype
from beartype.typing import Optional
from bokeh.io import show
from bokeh.layouts import column, row
from bokeh.models import (
    Button,
    ColumnDataSource,
    CustomJS,
    HoverTool,
    LassoSelectTool,
    Select,
)
from bokeh.plotting import figure
from deeporigin.data_hub.dataframe import DataFrame
from deeporigin.exceptions import DeepOriginException


@beartype
def scatter(
    df: DataFrame,
    x: Optional[str] = None,
    y: Optional[str] = None,
    size: Optional[str] = None,
):
    """function to make a scatter plot from a Deep Origin dataframe, with support for interactivity

    Args:
        df (DataFrame): Deep Origin dataframe
        x (Optional[str], optional): name of column to use for x axis. Defaults to None.
        y (Optional[str], optional): name of column to use for y axis. Defaults to None.
        size (Optional[str], optional): name of column to use for size. Defaults to None.

    `df` should be a Deep Origin dataframe with at least two numeric columns.

    """

    figure_width = 500
    select_width = int(figure_width * 0.3)
    js_code = _read_js_code()

    cols = df.attrs["metadata"]["cols"]
    numeric_cols = [col["name"] for col in cols if col["type"] in ["float", "integer"]]

    label_data = None

    select_cols = [
        col for col in df.attrs["metadata"]["cols"] if col["type"] == "select"
    ]
    if len(select_cols) > 0:
        select_col_names = [col["name"] for col in select_cols]
        select_col_options = [col["configSelect"]["options"] for col in select_cols]
        label_data = dict(zip(select_col_names, select_col_options))

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
    )

    # CDS for scatter plot
    scatter_source = ColumnDataSource(data)

    # CDS for marker
    marker_source = ColumnDataSource(_first_element_in_dict(data))

    # CDS to store lasso selection data
    lasso_selection_source = ColumnDataSource(dict(ids=[]))

    # CDS to store labels
    if label_data:
        label_source = ColumnDataSource(label_data)

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
        color="blue",
        alpha=0.6,
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
    if label_data:
        first_col = list(label_data.keys())[0]
        label_column_select = Select(
            title="Label Column",
            value=first_col,
            options=list(label_data.keys()),
            width=select_width,
        )

        label_select = Select(
            title="Label",
            value=label_data[first_col][0],
            options=list(set(label_data[first_col])),
            width=select_width,
        )

        # JavaScript callback to update the second select tool based on the selected column
        label_select_callback = CustomJS(
            args=dict(
                label_source=label_source,
                label_column_select=label_column_select,
                label_select=label_select,
            ),
            code="""
            // Get selected column name from first select widget
            const column = label_column_select.value;
            
            // Update options of the second select based on unique values of the selected column
            const column_data = label_source.data[column];
            const unique_values = Array.from(new Set(column_data));
            
            label_select.options = unique_values;
            label_select.value = unique_values[0];  // Set default to the first unique value
        """,
        )

        label_column_select.js_on_change("value", label_select_callback)

    # JavaScript callback to update data, axis labels, and point sizes on select change

    axes_callback = CustomJS(
        args=dict(
            scatter_source=scatter_source,
            x_select=x_select,
            y_select=y_select,
            size_select=size_select,
            df=df.to_dict("list"),
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

    # Button to access selected data from selected_source
    if label_data:
        label_button = Button(label="+Label")
        button_callback = CustomJS(
            args=dict(
                lasso_selection_source=lasso_selection_source,
                label_select=label_select,
                label_column_select=label_column_select,
            ),
            code=js_code["button_callback"],
        )
        label_button.js_on_click(button_callback)

    # Layout widgets and plot
    if label_data:
        layout = column(
            row(x_select, y_select, size_select),
            p,
            row(
                label_column_select,
                label_select,
                label_button,
            ),
        )
    else:
        layout = column(
            row(x_select, y_select, size_select),
            p,
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
