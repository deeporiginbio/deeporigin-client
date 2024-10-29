// Callback when a user presses the button to update labels in the scatter plot Bokeh figure

const selectedData = lasso_selection_source.data;
const label = label_select.value;
const labelColumn = label_column_select.value;

console.log("Will write to these rows:", selectedData.ids);
console.log("Will write this label:", label);
console.log("Will write to this column:", labelColumn);

const rowIds = selectedData.ids;

if (typeof window.deeporigin !== "undefined") {

    // Write the new value
    const updateChanges = rowIds.map(rowId => ({
        rowId: rowId,
        fieldChangeEvents: [{
            columnId: labelColumn,
            newValue: { selectedOptions: [label] } 
        }]
    }));

    deeporigin.dataHub.primaryDatabase.editRows({ changes: updateChanges });
}