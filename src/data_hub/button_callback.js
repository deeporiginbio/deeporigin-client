// Callback to update labels in the scatter plot Bokeh figure

const selectedData = lasso_selection_source.data;
const label = label_select.value;

console.log(scatter_source.data);

// Update colors for points in the lasso selection
scatter_source.data.id.forEach((id, i) => {
    if (selectedData.ids.includes(id)) {
        scatter_source.data.colors[i] = color_map[label];
        scatter_source.data.legend_labels[i] = label;
    }
});
scatter_source.change.emit();

const rowIds = selectedData.ids;

// Update external data if `deeporigin` is defined
if (window.deeporigin) {
    const updateChanges = rowIds.map(rowId => ({
        rowId,
        fieldChangeEvents: [
            {
                columnId: labelColumn,
                newValue: { selectedOptions: [label] }
            }
        ]
    }));

    deeporigin.dataHub.primaryDatabase.editRows({ changes: updateChanges });
}