// callback when a user hovers over a point in the scatter plot

if (cb_data.index.indices.length > 0) {
    // user is actually hovering over a point. 
    // bokeh has an annoying "feature" where you can hover
    // over more than 1 point, which we don't want.
    // we're going to arbitarily choose the first point
    var chosen_index = (cb_data.index.indices[0]);

    // update the red circle marker on all linked plots
    // to indicate the currently hovered point
    for (const key in marker_source.data) {
        marker_source.data[key][0] = scatter_source.data[key][chosen_index];
    }   

    marker_source.change.emit();

    // callback to update selection in database
    var id = scatter_source.data['id'][chosen_index];
    if (typeof window.deeporigin !== "undefined") {
        deeporigin.dataHub.primaryDatabase.clearRangeSelection();
        deeporigin.dataHub.primaryDatabase.addSelection({ selections: [{ rowId: id }] });
    
    }

};