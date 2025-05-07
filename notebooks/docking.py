

import marimo

__generated_with = "0.13.2"
app = marimo.App(width="medium")

with app.setup:
    import marimo as mo
    from deeporigin.drug_discovery import Protein, Ligand


@app.cell(hide_code=True)
def _():
    mo.md(
        """

    /// attention | Work in progress

    Functionality shown here is under active development. 
    ///

    
    # Docking using Deep Origin

    In this notebook we demonstrate how we can work with Proteins and Ligands to dock a Ligand to a Protein using the Deep Origin API.

        First, we construct a Protein object by specifying a PDB ID:
        """
    )
    return


@app.cell
def _():
    protein = Protein.from_pdb_id("1EBY")
    protein.show()
    return (protein,)


@app.cell(hide_code=True)
def _():
    mo.md(
    """

    ## Basic Protein Prep

    We notice that there are several heteratoms in the protein. We can remove them using the `remove_hetatm` method. This removes the waters and the ligand that is part of the structure we obtained from PDB:
    """
    
    )
    return


@app.cell
def _(protein):
    protein.remove_hetatm()
    protein.show()
    return


@app.cell(hide_code=True)
def _():
    mo.md(
    """
    ## Finding pockets

    We can now find pockets using the Pocket Finder tool"""
    
    )
    return


@app.cell
def _(protein):
    pockets = protein.find_pockets(pocket_count=1)
    pocket = pockets[0]
    return pocket, pockets


@app.cell(hide_code=True)
def _():
    mo.md(text=
         """
         We can inspect properties of the pocket by simply calling the pocket object:""")
    return


@app.cell
def _(pocket):
    pocket
    return


@app.cell(hide_code=True)
def _():
    mo.md("""
    ## Visualizing pocket

    We can visualize the pocket in the protein using the `show` method of the Protein class:""")
    return


@app.cell
def _(pockets, protein):
    protein.show(pockets=pockets)
    return


@app.cell(hide_code=True)
def _():
    mo.md(
        """
    ## Ligand to dock

    Here we construct the Ligand object to dock to this protein from a SMILES string:"""
    )
    return


@app.cell
def _():
    ligand = Ligand.from_smiles("CN(C)C(=O)c1cccc(-c2cn(C)c(=O)c3[nH]ccc23)c1")
    ligand.show()
    return (ligand,)


@app.cell(hide_code=True)
def _():
    mo.md(

    """

    ## Docking

    We can now dock the ligand to the protein using the `dock` method of the Protein class. This will return a SDF file with the poses of the ligand in the protein. We can visualize these poses using the `show` method of the Protein class:

    """
    )
    return


@app.cell
def _(ligand, pocket, protein):
    poses_sdf = protein.dock(ligand=ligand,pocket=pocket)
    protein.show(sdf_file=poses_sdf)
    return


@app.cell(hide_code=True)
def _():
    return


if __name__ == "__main__":
    app.run()
