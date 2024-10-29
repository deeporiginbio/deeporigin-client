// callback on lasso selection

function debounce(func, wait) {
    let timeout;
    return function (...args) {
        const context = this;
        clearTimeout(timeout);
        timeout = setTimeout(() => func.apply(context, args), wait);
    };
}

const processSelection = () => {
    const selected_indices = scatter_source.selected.indices;
    let selected_data = selected_indices.map(i => {
        let row = {};
        for (const key in scatter_source.data) {
            row[key] = scatter_source.data[key][i];
        }
        return row;
    });

    const ids = selected_data.map(item => item.id);
    console.log("Lasso selected points data:", ids);

    lasso_selection_source.data.ids = ids;
    lasso_selection_source.change.emit();
    
};


const debouncedSelection = debounce(processSelection, 50);
debouncedSelection();