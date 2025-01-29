"""this module contains high level utility functions that makes it easier to use tools by providing wrappers over tools"""

from typing import Literal

from beartype import beartype
from deeporigin.data_hub import api
from deeporigin.exceptions import DeepOriginException
from deeporigin.tools import run
from deeporigin.utils.config import construct_resource_url

VINA_DB = "autodock-vina"
LIGAND_PREP_DB = "ligand-prep-meeko"
ABFE_DB = "ABFE"


charge_methods = Literal["gas", "bcc"]
lig_force_fields = Literal["gaff2", "openff"]


@beartype
def _ensure_database(name: str) -> dict:
    """ensure that a database exists with the given name. If it doesn't exist, create it"""

    databases = api.list_rows(row_type="database")

    database = [db for db in databases if db["hid"] == name]

    if len(database) == 0:
        # make a new DB
        print("🧬 Creating a database called {name}...")
        database = api.create_database(
            hid=name,
            hid_prefix=name,
            name=name,
        )

    database = api.describe_database(database_id=name)
    return database


@beartype
def _ensure_columns(
    *,
    database: dict,
    required_columns: list[dict],
):
    """ensure that columns exist with the given names (and types). If they don't exist, create them"""

    existing_column_names = []
    if "cols" in list(database.keys()):
        existing_column_names = [col["name"] for col in database.cols]

    # check if we need to make columns
    for item in required_columns:
        column_name = item["name"]
        column_type = item["type"]

        if column_name in existing_column_names:
            continue
        print(f"🧬 Making column named: {column_name} in {database.hid}")

        api.add_database_column(
            cardinality="one",
            database_id=database.hid,
            name=column_name,
            required=False,
            type=column_type,
        )


@beartype
def _ensure_db_for_abfe() -> dict:
    """ensure that there is a database for FEP on Data Hub"""

    required_columns = [
        dict(name="ligand", type="file"),
        dict(name="protein", type="file"),
        dict(name="prep_output", type="file"),
        dict(name="emeq_output", type="file"),
    ]

    database = _ensure_database(ABFE_DB)

    _ensure_columns(
        database=database,
        required_columns=required_columns,
    )

    return database


def _abfe_system_prep(
    *,
    ligand_file: str,
    protein_file: str,
    keep_waters: bool = True,
    save_gmx_files: bool = False,
    is_lig_protonated: bool = True,
    is_protein_protonated: bool = True,
    do_loop_modelling: bool = False,
    charge_method: charge_methods = "bcc",
    lig_force_field: lig_force_fields = "gaff2",
    padding: float = 1.0,
    system_name: str = "complex",
):
    """high level function to prepare ligand and protein files for FEP on Data Hub"""

    # input validation
    if padding < 0.5 or padding > 2:
        raise DeepOriginException("Padding must be greater than 0.5, and less than 2")

    database = _ensure_db_for_abfe()

    url = construct_resource_url(
        name=ABFE_DB,
        row_type="database",
    )

    print(f"Using database at: {url}")

    print(f"🧬 Print uploading files to database called: {database.hid}...")

    response = api.upload_file_to_new_database_row(
        database_id=ABFE_DB,
        column_id="ligand",
        file_path=ligand_file,
    )

    row_id = response.rows[0].hid

    file = api.upload_file(file_path=protein_file)

    api.assign_files_to_cell(
        file_ids=[file.id],
        database_id=ABFE_DB,
        column_id="protein",
        row_id=row_id,
    )

    tool_id = "deeporigin/md-suite-prep"
    inputs = {
        "ligand": {
            "columnId": "ligand",
            "rowId": row_id,
            "databaseId": database.hid,
        },
        "protein": {
            "columnId": "protein",
            "rowId": row_id,
            "databaseId": database.hid,
        },
        "force": 1,
        "test_run": 0,
        "system": system_name,
        "include_ligands": 1,
        "include_protein": 1,
        "sysprep_params": {
            "is_protein_protonated": is_protein_protonated,
            "do_loop_modelling": do_loop_modelling,
            "force_field": "ff14SB",
            "padding": padding,
            "keep_waters": keep_waters,
            "save_gmx_files": save_gmx_files,
            "is_lig_protonated": is_lig_protonated,
            "charge_method": charge_method,
            "lig_force_field": lig_force_field,
        },
    }

    outputs = {
        "output_file": {
            "columnId": "prep_output",
            "rowId": row_id,
            "databaseId": database.hid,
        }
    }

    run._process_job(
        inputs=inputs,
        outputs=outputs,
        tool_id=tool_id,
        cols=database.cols,
    )


@beartype
def _ensure_db_for_vina_outputs() -> None:
    """makes a DB for outputs of vina"""

    databases = api.list_rows(row_type="database")

    cols = ["docked-ligand", "scores", "SMILES"]
    col_types = ["file", "file", "text"]

    database = [db for db in databases if db["hid"] == VINA_DB]

    if len(database) == 0:
        # make a new DB
        database = api.create_database(
            hid=VINA_DB,
            hid_prefix="adv",
            name=VINA_DB,
        )

    # get info about database
    database = api.describe_database(database_id=VINA_DB)

    existing_cols = []
    if "cols" in list(database.keys()):
        existing_cols = [col["name"] for col in database.cols]

    # check if we need to make columns

    for col, col_type in zip(cols, col_types):
        if col in existing_cols:
            continue
        print(f"making {col}")

        api.add_database_column(
            cardinality="one",
            database_id=VINA_DB,
            name=col,
            required=False,
            type=col_type,
        )


@beartype
def _ensure_db_for_ligand_prep_meeko() -> None:
    databases = api.list_rows(row_type="database")

    cols = ["input-ligand", "output-ligand"]

    database = [db for db in databases if db["hid"] == LIGAND_PREP_DB]

    if len(database) == 0:
        # make a new DB
        database = api.create_database(
            hid=LIGAND_PREP_DB,
            hid_prefix="lpm",
            name=LIGAND_PREP_DB,
        )

    # get info about database
    database = api.describe_database(database_id=LIGAND_PREP_DB)

    existing_cols = []
    if "cols" in list(database.keys()):
        existing_cols = [col["name"] for col in database.cols]

    # check if we need to make columns
    for col in cols:
        if col in existing_cols:
            continue
        print(f"making {col}")

        api.add_database_column(
            cardinality="one",
            database_id=LIGAND_PREP_DB,
            name=col,
            required=False,
            type="file",
        )


@beartype
def smiles_to_sdf(smiles: str, sdf_path: str) -> None:
    """convert a SMILES string to a SDF file"""

    from rdkit import Chem
    from rdkit.Chem import AllChem, SDWriter

    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        print(f"Invalid SMILES: {smiles}")

    try:
        Chem.Kekulize(mol)
    except ValueError:
        print(f"Failed to kekulize: {smiles}")

    mol = Chem.AddHs(mol)

    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    AllChem.UFFOptimizeMolecule(mol)

    with SDWriter(sdf_path) as writer:
        writer.write(mol)


@beartype
def csv_to_sdfs(
    csv_file: str,
    smiles_column_name: str,
) -> list:
    """given a CSV file, convert all SMILES strings in that CSV file to SDF files"""

    import pandas as pd

    df = pd.read_csv(csv_file)
    smiles_list = list(df[smiles_column_name])

    all_sdf_paths = []

    for i, smiles in enumerate(smiles_list):
        sdf_path = f"{str(i)}.sdf"
        smiles_to_sdf(
            smiles,
            sdf_path=sdf_path,
        )
        all_sdf_paths.append(sdf_path)
    return all_sdf_paths


@beartype
def prep_ligands(*, csv_file: str, column_name: str):
    """high level function to prepare ligands from a CSV file with SMILES strings using Lingand Prep Meeko Tool on Deep Origin"""

    _ensure_db_for_ligand_prep_meeko()

    sdf_files = csv_to_sdfs(csv_file, column_name)

    row_ids = []

    for sdf_file in sdf_files:
        print(f"Uploading {sdf_file}...")
        response = api.upload_file_to_new_database_row(
            database_id=LIGAND_PREP_DB,
            column_id="input-ligand",
            file_path=sdf_file,
        )

        row_ids.append(response.rows[0].hid)

    run_ligand_prep_on_db(
        database_id=LIGAND_PREP_DB,
        input_column_name="input-ligand",
        output_column_name="output-ligand",
    )


def run_ligand_prep_on_db(
    *,
    database_id: str,
    input_column_name: str,
    output_column_name: str,
):
    """utility function to run Ligand Prep (Meeko) on all ligands in a database"""

    database = api.describe_database(database_id=database_id)
    rows = api.list_database_rows(database_row_id=database_id)
    col_mapper = {col["name"]: col["id"] for col in database.cols}

    for row in rows:
        if "fields" in row.keys():
            columns_with_files = [field["columnId"] for field in row.fields]
        else:
            columns_with_files = []

        if (
            col_mapper["input-ligand"] in columns_with_files
            and col_mapper["output-ligand"] not in columns_with_files
        ):
            run.ligand_prep(
                database_id=database_id,
                ligand_column_name=input_column_name,
                output_column_name=output_column_name,
                row_id=row.hid,
            )


@beartype
def dock_ligands_to_protein(
    *,
    ligand_database_name: str,
    ligand_column_name: str,
    protein_database_name: str,
    protein_column_name: str,
    protein_row: str,
    search_space: dict,
    docking: dict,
):
    """Dock all ligands in a ligand database to a prepped protein in a protein database, and store those results in a new database"""

    _ensure_db_for_vina_outputs()

    # figure out which rows of the ligand DB are valid
    ligand_db = api.describe_database(database_id=ligand_database_name)
    protein_db = api.describe_database(database_id=protein_database_name)
    vina_db = api.describe_database(database_id=VINA_DB)
    ligand_rows = api.list_database_rows(database_row_id=ligand_database_name)

    if len(ligand_rows) == 0:
        raise ValueError("No rows in ligand DB")

    col_mapper = {col["name"]: col["id"] for col in ligand_db.cols}

    rows_with_ligand = []

    for row in ligand_rows:
        columns_with_files = [field["columnId"] for field in row.fields]

        if col_mapper[ligand_column_name] in columns_with_files:
            rows_with_ligand.append(row.hid)

    # make new rows, as many as there are rows_with_ligand
    n_rows = len(rows_with_ligand)
    if n_rows == 0:
        raise ValueError(
            "There are no valid ligands that can be used as inputs to Vina in the ligand DB"
        )

    print(f"Making {n_rows} new rows in the {VINA_DB} database to hold outputs...")
    response = api.make_database_rows(
        database_id=VINA_DB,
        n_rows=n_rows,
    )

    output_rows = [row["hid"] for row in response.rows]

    # for each row with ligand, run a vina job
    for ligand_row, output_row in zip(
        rows_with_ligand,
        output_rows,
    ):
        inputs = dict(
            receptor={
                "rowId": protein_row,
                "columnId": protein_column_name,
                "databaseId": protein_database_name,
            },
            ligand={
                "rowId": ligand_row,
                "columnId": ligand_column_name,
                "databaseId": ligand_database_name,
            },
            searchSpace=search_space,
            docking=docking,
        )
        outputs = dict(
            output_file={
                "rowId": output_row,
                "columnId": "docked-ligand",
                "databaseId": VINA_DB,
            },
            scores_file={
                "rowId": output_row,
                "columnId": "scores",
                "databaseId": VINA_DB,
            },
        )

        run._process_job(
            cols=protein_db.cols + ligand_db.cols + vina_db.cols,
            inputs=inputs,
            outputs=outputs,
            tool_id="deeporigin/autodock-vina",
        )
