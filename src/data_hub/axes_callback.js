// JS callback to handle axes selection
// this is called every time the user makes a selection
// on what column of the dataframe to plot onto what axes


const x_data = df[x_select.value];
const y_data = df[y_select.value];
const size_data = df[size_select.value];


// Normalize size data for better visualization (e.g., map values to a range)
const min_size = Math.min(...size_data);
const max_size = Math.max(...size_data);

// Define a scaling function for sizes
const sizes = size_data.map(value => {
    return 2 + 15 * (value - min_size) / (max_size - min_size);
});

// Update the data source
scatter_source.data.x = x_data;
scatter_source.data.y = y_data;
scatter_source.data.size = sizes;

// Update the axis labels
x_axis.axis_label = x_select.value;
y_axis.axis_label = y_select.value;

// Trigger data change
scatter_source.change.emit();
