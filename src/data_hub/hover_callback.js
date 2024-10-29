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