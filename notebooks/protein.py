

import marimo

__generated_with = "0.13.2"
app = marimo.App(width="medium")

with app.setup:
    import marimo as mo


@app.cell(hide_code=True)
def _():
    mo.md(
        """
        /// attention | Work in progress

        Functionality shown here is under active development. 
        ///

        # Working with proteins

        This notebook shows you how to use the Protein class in the Deep Origin Drug Discovery toolbox
        """
    )
    return


@app.cell
def _():
    from deeporigin.drug_discovery import Complex, EXAMPLE_DATA_DIR, Protein, Ligand
    return Complex, EXAMPLE_DATA_DIR, Protein


@app.cell(hide_code=True)
def _():
    mo.md(
        """
        ## Making a Protein

        We can construct Proteins using a number of ways. Proteins can be pulled directly from the PDB by specifying a PDB ID. 

        Here, we pull a conotoxin. Displaying the object shows a helpful widget with information about the protein
        """
    )
    return


@app.cell
def _(Protein):
    conotoxin = Protein.from_pdb_id("2JUQ")
    conotoxin
    return (conotoxin,)


@app.cell(hide_code=True)
def _():
    mo.md(r"""We can use the `show` method of the Protein class to visualize its 3D structure""")
    return


@app.cell
def _(conotoxin):
    conotoxin.show()
    return


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""
        ## Making a protein from a file

        Now we demonstrate how we can construct a Protein object from a PDB file that is stored locally
        """
    )
    return


@app.cell
def _(Complex, EXAMPLE_DATA_DIR):
    sim = Complex.from_dir(EXAMPLE_DATA_DIR)
    return (sim,)


@app.cell
def _(sim):
    protein = sim.protein
    protein.show()
    return (protein,)


@app.cell(hide_code=True)
def _():
    mo.md(
        r"""
        ## Operations on proteins

        The `Protein` class has several operations. For example, we can remove all waters from the protein
        """
    )
    return


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


if __name__ == "__main__":
    app.run()
