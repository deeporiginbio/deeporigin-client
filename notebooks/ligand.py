

import marimo

__generated_with = "0.13.2"
app = marimo.App(width="medium")


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
        # Working with ligands

        This notebook shows you how to use the Ligand class in the Deep Origin Drug Discovery toolbox
        """
    )
    return


@app.cell
def _():
    from deeporigin.drug_discovery import EXAMPLE_DATA_DIR, Ligand
    return EXAMPLE_DATA_DIR, Ligand


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
        ## Making a Ligand

        ### From an identifier

        We can make Ligands using an identifier. Molecules are pulled from PubChem using [pubchempy](https://github.com/mcs07/PubChemPy)
        """
    )
    return


@app.cell
def _(Ligand):
    ligand = Ligand.from_identifier("ATP")
    ligand.show()
    return (ligand,)


@app.cell(hide_code=True)
def _(mo):
    mo.md(r"""### from SMILES Strings""")
    return


@app.cell
def _(Ligand):
    pcm = Ligand.from_smiles("CC(=O)Nc1ccc(O)cc1")
    pcm.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        r"""
        ### From a file

        We can make Ligands from a SDF file
        """
    )
    return


@app.cell
def _(EXAMPLE_DATA_DIR, Ligand):
    brd2 = Ligand.from_sdf(EXAMPLE_DATA_DIR / "brd-2.sdf")
    brd2.show()
    return


@app.cell(hide_code=True)
def _(mo):
    mo.md(
        """
        ## Ligand operations

        The Ligand class offers many functions. For example, we can use the DO Molecular property predictor to predict ADMET properties:
        """
    )
    return


@app.cell
def _(ligand):
    ligand.admet_properties()
    return


app._unparsable_cell(
    r"""
    We can inspect the ligand's properties, such as SMILES string:


    """,
    column=None, disabled=False, hide_code=True, name="_"
)


@app.cell
def _(ligand):
    ligand.properties
    return


@app.cell
def _():
    import marimo as mo
    return (mo,)


if __name__ == "__main__":
    app.run()
