

import marimo

__generated_with = "0.13.2"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Jobs

        This notebook shows you how to use the Job dashboard to query jobs and make plots of various things
        """
    )
    return


@app.cell
def _():
    from bokeh.io import output_notebook, show
    from bokeh.models import ColumnDataSource, Jitter
    from bokeh.palettes import Category10, Category20
    from bokeh.plotting import figure
    from bokeh.transform import factor_cmap, jitter
    import pandas as pd

    from deeporigin.tools.job import get_dataframe

    output_notebook()
    return (
        Category10,
        Category20,
        ColumnDataSource,
        factor_cmap,
        figure,
        get_dataframe,
        jitter,
        show,
    )


@app.cell
def _(get_dataframe):
    def get_jobs(tool_key="deeporigin.abfe-end-to-end"):
        df = get_dataframe(tool_key=tool_key, only_with_status=["Succeeded"])
        df = df[(df["test_run"] == False)
        ]
        return df
    return (get_jobs,)


@app.cell
def _(get_jobs):
    get_jobs()
    return


@app.cell
def _():
    from deeporigin.data_hub import api

    df_proteins = api.get_dataframe("Proteins")
    df_proteins = df_proteins.drop(columns=["Validation Status"])
    mapper = df_proteins.to_dict()["Protein"]
    # Transform the mapper dictionary:
    # 1. Remove ".pdb" from string values.
    # 2. Remove "protein" from string values.
    # 3. Split string values at the first underscore and keep the first part.
    mapper = {
        k: (
            v.replace(".pdb", "").replace("protein", "").split("_", 1)[0]
            if isinstance(v, str) else v
        )
        for k, v in mapper.items()
    }
    return (mapper,)


@app.cell
def _(
    Category10,
    Category20,
    ColumnDataSource,
    factor_cmap,
    figure,
    get_jobs,
    jitter,
    mapper,
    show,
):
    def make_plot(color_by: str, group_by: str):

        df = get_jobs()
    
        df[group_by] = (
            df[group_by]
            .map(mapper)  # apply your mapping
            .fillna(df[group_by])  # fallback to original id
            .fillna("None")  # handle NaNs
        )

        df[color_by] = (
            df[color_by]
            .map(mapper)  # apply your mapping
            .fillna(df[color_by])  # fallback to original id
            .fillna("None")  # handle NaNs
        )
    
        source = ColumnDataSource(df)
    
        proteins = list(df[group_by].unique())
        tool_versions = sorted(df[color_by].unique())  # for coloring
    
        n_versions = len(tool_versions)
        palette = (
            Category10[min(10, n_versions)]
            if n_versions <= 10
            else Category20[min(20, n_versions)]
        )
    
        p = figure(
            width=700,
            height=400,
            title="Run duration vs. protein (colored by tool version)",
            x_axis_label="Protein",
            y_axis_label="Run duration (minutes)",
            x_range=proteins,
            y_range=(0, 1500),  # keep your previous y-range
            tools="pan,wheel_zoom,box_zoom,reset,save",
            toolbar_location="above",
        )
    
        p.scatter(
            x=jitter(group_by, width=0.2, range=p.x_range),
            y="run_duration_minutes",
            source=source,
            size=8,
            marker="circle",
            fill_alpha=0.7,
            line_color=None,
            color=factor_cmap(color_by, palette=palette, factors=tool_versions),
            legend_field=color_by,  # legend now shows tool versions
        )
    
        p.legend.title = color_by
        p.legend.location = "bottom_left"
    
        show(p)


    return (make_plot,)


@app.cell
def _(make_plot):
    make_plot(color_by="tool_version", group_by="protein_id")
    return


@app.cell
def _(make_plot):
    make_plot(color_by="protein_id", group_by="tool_version")
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


@app.cell
def _():
    return


if __name__ == "__main__":
    app.run()
