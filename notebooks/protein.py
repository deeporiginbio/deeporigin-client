

import marimo

__generated_with = "0.13.2"
app = marimo.App(width="medium")


app._unparsable_cell(
    r"""
    # Protein

    This notebook shows you how to use the Protein class in the Deep Origin Drug Discovery toolbox
    """,
    column=None, disabled=False, hide_code=True, name="_"
)


@app.cell
def _():
    from deeporigin.drug_discovery import Complex, EXAMPLE_DATA_DIR, Protein, Ligand
    return Complex, EXAMPLE_DATA_DIR


@app.cell
def _(Complex, EXAMPLE_DATA_DIR):
    sim = Complex.from_dir(EXAMPLE_DATA_DIR)
    return (sim,)


@app.cell
def _(sim):
    protein = sim.protein
    protein.show()
    return (protein,)


@app.cell
def _(sim):
    sim.protein.remove_water()
    sim.protein.show()
    return


@app.cell
def _(protein):
    # Remove water molecules (HOH) and ions (NA, CL)
    protein.remove_resnames(exclude_resnames=["HOH", "NA", "CL"])
    return


@app.cell
def _(protein):
    protein.show()
    return


@app.cell
def _():
    import marimo as mo
    return


if __name__ == "__main__":
    app.run()
