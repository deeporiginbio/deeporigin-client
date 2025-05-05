

import marimo

__generated_with = "0.13.2"
app = marimo.App(width="medium")


app._unparsable_cell(
    r"""
    # ABFE on Deep Origin

    This notebook demonstrates how one runs ABFE on Deep Origin
    """,
    name="_"
)


@app.cell
def _():
    from deeporigin.drug_discovery import Complex, EXAMPLE_DATA_DIR, Protein, Ligand
    return Complex, EXAMPLE_DATA_DIR


@app.cell
def _(Complex, EXAMPLE_DATA_DIR):
    sim = Complex.from_dir(EXAMPLE_DATA_DIR)
    sim.connect()
    return (sim,)


@app.cell
def _(sim):
    sim.abfe.show_results()
    return


@app.cell
def _():
    import marimo as mo
    return


if __name__ == "__main__":
    app.run()
