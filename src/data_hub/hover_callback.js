// Callback when a user hovers over a point in the scatter plot
if (cb_data.index && cb_data.index.indices && cb_data.index.indices.length > 0) {
    // User is actually hovering over a point. 
    // Bokeh allows hovering over multiple points; we select the first.
    const chosenIndex = cb_data.index.indices[0];

    // Update the red circle marker on all linked plots to indicate the hovered point
    for (const key in marker_source.data) {
        if (scatter_source.data[key]) {
            marker_source.data[key][0] = scatter_source.data[key][chosenIndex];
        }
    }

    marker_source.change.emit();


    // now update the table
    const hidden_cols = ["Validation Status", "colors", "legend_labels"];
    const data = scatter_source.data;
    const columns = Object.keys(data);
    const columnNames = [];
    const values = [];
    for (const col of columns) {
        if (hidden_cols.includes(col)) { continue; }

        columnNames.push(col);
        // round floats so they don't look ugly
        var value = data[col][chosenIndex];

        if (typeof value === 'number' && !Number.isInteger(value)) {
            value = value.toFixed(3);
        }
        values.push(value);
    }
    table_source.data = { 'Column Name': columnNames, 'Value': values };
    table_source.change.emit();

    // Callback to update selection in the database
    const id = scatter_source.data['id'] ? scatter_source.data['id'][chosenIndex] : null;
    if (id && typeof window.deeporigin !== "undefined") {
        try {
            deeporigin.dataHub.primaryDatabase.clearRangeSelection();
            deeporigin.dataHub.primaryDatabase.addSelection({ selections: [{ rowId: id }] });
        } catch (error) {
            console.error("Error updating selection in database:", error);
        }
    }

    console.log("Hovered point ID:", id);
}