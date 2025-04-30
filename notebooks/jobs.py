

import marimo

__generated_with = "0.13.2"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ## Jobs

        This notebook shows you how to use the Job dashboard
        """
    )
    return


@app.cell
def _():
    from deeporigin.tools.job import get_dataframe
    return (get_dataframe,)


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
    df = get_jobs("deeporigin.bulk-docking")
    df
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
