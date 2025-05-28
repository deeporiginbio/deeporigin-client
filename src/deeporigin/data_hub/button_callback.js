// Callback to update labels in the scatter plot Bokeh figure
// when the button is pressed

// get the points within the lasso selection
const selectedData = lasso_selection_source.data;

// get the label we need to assign
const label = label_select.value;

console.log("Assigning label:", label);

// Update colors for points in the lasso selection
scatter_source.data.id.forEach((id, i) => {
    if (selectedData.ids.includes(id)) {
        scatter_source.data.colors[i] = color_map[label];
        scatter_source.data.legend_labels[i] = label;
    }
});
scatter_source.change.emit();

const rowIds = selectedData.ids;

console.log("Assigning to rows:", rowIds);

// write changes to the deep origin database
if (window.deeporigin) {
    const updateChanges = rowIds.map(rowId => ({
        rowId,
        fieldChangeEvents: [
            {
                columnId: label_column,
                newValue: { selectedOptions: [label] }
            }
        ]
    }));

    deeporigin.dataHub.primaryDatabase.editRows({
        checkPreviousValue: false, 
        changes: updateChanges
    });
}