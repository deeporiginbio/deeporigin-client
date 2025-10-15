import marimo

__generated_with = "0.16.5"
app = marimo.App(width="medium")

with app.setup:
    import marimo as mo

    from deeporigin.drug_discovery import BRD_DATA_DIR, Complex


@app.cell(hide_code=True)
def _():
    mo.md(
        """
    /// attention | Work in progress

    Functionality shown here is under active development. 
    ///


    # ABFE on Deep Origin

    This notebook demonstrates how one runs ABFE on Deep Origin

    ## Initialization 

    We first create a Complex from a protein and a set of ligands
    """
    )
    return


@app.cell
def _():
    sim = Complex.from_dir(BRD_DATA_DIR)
    return (sim,)


@app.cell
def _(sim):
    sim.show_ligands()
    return


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""
    ## System preparation

    Before we start a ABFE run, we can prepare the system using `sim.prepare`. Doing so shows a preview of the prepared system.
    """
    )
    return


@app.cell
def _(sim):
    sim.prepare(ligand=sim.ligands[0])
    return


@app.cell(hide_code=True)
def _():
    mo.md(
        """
    ## Job control

    In this section we show how we can

    - start a job
    - view status of a job
    - cancel a job


    ### Starting a ABFE run
    """
    )
    return


@app.cell
def _(sim):
    sim.abfe.set_test_run()
    jobs = sim.abfe.run_end_to_end(ligands=[sim.ligands[0]])
    this_job = jobs[0]
    return (this_job,)


@app.cell(hide_code=True)
def _():
    mo.md(
        """
    ### Inspecting a single run

    In this section we show how we can use the Job class to inspect a the run we just started:
    """
    )
    return


@app.cell
def _(this_job):
    this_job
    return


@app.cell(hide_code=True)
def _():
    mo.md(
        """
    ### Cancelling a job

    In this section we cancel the job we just started
    """
    )
    return


@app.cell
def _(this_job):
    this_job.cancel()
    return


@app.cell(hide_code=True)
def _():
    mo.md(
        """
    ## Viewing results

    Here we view results from previous (completed) runs
    """
    )
    return


@app.cell
def _(sim):
    sim.abfe.show_results()
    return


@app.cell(hide_code=True)
def _():
    mo.md(
        """
    ## Viewing all jobs

    In this section we dig deeper into viewing jobs.
    """
    )
    return


@app.cell
def _(sim):
    df = sim.abfe.show_jobs(summary=True)
    df
    return


@app.cell(hide_code=True)
def _():
    mo.md(
        """
    ## Viewing all jobs in an org

    We can view all jobs in an org, so we can filter post-hoc by user, tool, status, etc
    """
    )
    return


@app.cell
def _():
    from deeporigin.tools.job import get_dataframe

    get_dataframe()
    return


@app.cell
def _():
    from deeporigin.platform import Client

    client = Client()
    return (client,)


@app.cell
def _(client, sim):
    sim2 = Complex(protein=sim.protein, client=client)
    return (sim2,)


@app.cell
def _(sim2):
    sim2.client
    return


if __name__ == "__main__":
    app.run()
