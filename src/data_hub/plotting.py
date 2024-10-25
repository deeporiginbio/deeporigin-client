from bokeh.io import show
from bokeh.layouts import column, row
from bokeh.models import ColumnDataSource, CustomJS, Select
from bokeh.plotting import figure
from deeporigin.data_hub.dataframe import DataFrame


def create_interactive_scatter(df: DataFrame):
    """function to make a scatter plot from a Deep Origin dataframe"""

    figure_width = 500

    if df.shape[1] < 2:
        raise ValueError("DataFrame must contain at least two columns.")

    cols = df.attrs["metadata"]["cols"]
    cols = [col["name"] for col in cols if col["type"] in ["float", "integer"]]

    # Set up initial data for scatter plot
    initial_x = cols[0]
    initial_y = cols[1]
    initial_size = cols[0]

    # Create a ColumnDataSource with initial data
    source = ColumnDataSource(
        data={
            "x": df[initial_x],
            "y": df[initial_y],
            "size": df[initial_size],
        }
    )

    # Create the scatter plot figure
    p = figure(
        width=figure_width,
        height=figure_width,
        title="Deep Origin / " + df.attrs["metadata"]["name"],
    )
    p.toolbar.logo = None

    p.scatter(
        x="x",
        y="y",
        source=source,
        size="size",
        color="blue",
        alpha=0.6,
    )

    # Set initial axis labels
    p.xaxis.axis_label = initial_x
    p.yaxis.axis_label = initial_y

    # Create dropdown selectors for X and Y axes
    x_select = Select(
        title="X-Axis",
        value=initial_x,
        options=cols,
        width=int(figure_width * 0.3),
    )
    y_select = Select(
        title="Y-Axis",
        value=initial_y,
        options=cols,
        width=int(figure_width * 0.3),
    )
    size_select = Select(
        title="Size",
        value=initial_size,
        options=cols,
        width=int(figure_width * 0.3),
    )

    # JavaScript callback to update data, axis labels, and point sizes on select change
    callback = CustomJS(
        args=dict(
            source=source,
            x_select=x_select,
            y_select=y_select,
            size_select=size_select,
            df=df.to_dict("list"),
            x_axis=p.xaxis[0],
            y_axis=p.yaxis[0],
        ),
        code="""
        const x_data = df[x_select.value];
        const y_data = df[y_select.value];
        const size_data = df[size_select.value];
        
        // Normalize size data for better visualization (e.g., map values to a range)
        const min_size = Math.min(...size_data);
        const max_size = Math.max(...size_data);
        
        // Define a scaling function for sizes
        const sizes = size_data.map(value => {
            return 2 + 15 * (value - min_size) / (max_size - min_size); // Scale sizes between 10 and 40
        });
        
        // Update the data source
        source.data = {'x': x_data, 'y': y_data, 'size': sizes};
        
        // Update the axis labels
        x_axis.axis_label = x_select.value;
        y_axis.axis_label = y_select.value;
        
        // Trigger data change
        source.change.emit();
    """,
    )

    # Attach the callback to the select widgets
    x_select.js_on_change("value", callback)
    y_select.js_on_change("value", callback)
    size_select.js_on_change("value", callback)

    # Layout widgets and plot
    layout = column(
        row(x_select, y_select, size_select),
        p,
    )
    show(layout)
