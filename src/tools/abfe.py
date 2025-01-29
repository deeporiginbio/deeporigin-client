from typing import Literal

from beartype import beartype
from deeporigin.data_hub import api
from deeporigin.exceptions import DeepOriginException
from deeporigin.tools import run
from deeporigin.tools.toolkit import _ensure_columns, _ensure_database
from deeporigin.utils.config import construct_resource_url

ABFE_DB = "ABFE"
charge_methods = Literal["gas", "bcc"]
lig_force_fields = Literal["gaff2", "openff"]
force_fields = Literal["ff14SB", "ff99SB-ildn"]


@beartype
def init(
    *,
    ligand_file: str,
    protein_file: str,
) -> str:
    """Initialize an ABFE run. Upload ligand and protein files to Data Hub."""

    _ensure_db_for_abfe()

    url = construct_resource_url(
        name=ABFE_DB,
        row_type="database",
    )

    print(f"Using database at: {url}")
    print("🧬 Uploading files to database...")

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

    print(f"🧬 Files uploaded to row {row_id}.")

    return row_id


@beartype
def _ensure_db_for_abfe() -> dict:
    """ensure that there is a database for FEP on Data Hub"""

    required_columns = [
        dict(name="ligand", type="file"),
        dict(name="protein", type="file"),
        dict(name="complex_prep_output", type="file"),
        dict(name="solvation_prep_output", type="file"),
        dict(name="emeq_output", type="file"),
    ]

    database = _ensure_database(ABFE_DB)

    _ensure_columns(
        database=database,
        required_columns=required_columns,
    )

    return database


@beartype
def emeq(
    *,
    row_id: str,
    system_name: str = "complex",
) -> None:
    tool_id = "deeporigin/md-suite-emeq"

    database = _ensure_database(ABFE_DB)

    inputs = {
        "input": {
            "columnId": "complex_prep_output",
            "rowId": row_id,
            "databaseId": database.hid,
        },
        "force": 1,
        "test_run": 0,
        "system": system_name,
        "run_name": "test-run",
        "threads": 0,
        "em_solvent": True,
        "em_all": True,
        "nvt_heating_ns": 0.1,
        "npt_reduce_restraints_ns": 0.2,
        "from_run": "__USE_SYSTEM",
        "emeq_md_options": {
            "Δt": 0.004,
            "T": 298.15,
            "cutoff": 0.9,
            "fourier_spacing": 0.12,
            "hydrogen_mass": 2,
        },
    }

    outputs = {
        "output_file": {
            "columnId": "emeq_output",
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
def complex_prep(
    *,
    row_id: str,
    keep_waters: bool = True,
    save_gmx_files: bool = False,
    is_lig_protonated: bool = True,
    is_protein_protonated: bool = True,
    do_loop_modelling: bool = False,
    charge_method: charge_methods = "bcc",
    lig_force_field: lig_force_fields = "gaff2",
    padding: float = 1.0,
    system_name: str = "complex",
    force_field: force_fields = "ff14SB",
) -> None:
    """Function to prepare uploaded Ligand and protein files using Deep Origin MDSuite. Use this function to run system prep on a ligand and protein pair, that exist as files on a row in the ABFE database."""
    kwargs = locals()

    kwargs["include_protein"] = 1
    kwargs["include_ligands"] = 1
    kwargs["output_column_name"] = "complex_prep_output"

    _prep(**kwargs)


@beartype
def solvation_prep(
    *,
    row_id: str,
    is_lig_protonated: bool = True,
    charge_method: charge_methods = "bcc",
    lig_force_field: lig_force_fields = "gaff2",
    system_name: str = "ligand_prep",
    force_field: force_fields = "ff14SB",
) -> None:
    """Function to prepare uploaded Ligand and protein files using Deep Origin MDSuite. Use this function to run system prep on a ligand and protein pair, that exist as files on a row in the ABFE database."""
    kwargs = locals()

    kwargs["include_protein"] = 0
    kwargs["include_ligands"] = 1
    kwargs["output_column_name"] = "solvation_prep_output"

    _prep(**kwargs)


@beartype
def _prep(
    *,
    row_id: str,
    include_protein: int,
    include_ligands: int,
    system_name: str,
    output_column_name: str,
    keep_waters: bool = True,
    save_gmx_files: bool = False,
    is_lig_protonated: bool = True,
    is_protein_protonated: bool = True,
    do_loop_modelling: bool = False,
    charge_method: charge_methods = "bcc",
    lig_force_field: lig_force_fields = "gaff2",
    padding: float = 1.0,
    force_field: force_fields = "ff14SB",
) -> None:
    """Function to prepare uploaded Ligand and protein files using Deep Origin MDSuite. Use this function to run system prep on a ligand and protein pair, that exist as files on a row in the ABFE database.

    Args:
        row_id (str): row id of the ligand and protein files.
        keep_waters (bool, optional): whether to keep waters. Defaults to True.
        save_gmx_files (bool, optional): whether to save gmx files. Defaults to False.
        is_lig_protonated (bool, optional): whether the ligand is protonated. Defaults to True.
        is_protein_protonated (bool, optional): whether the protein is protonated. Defaults to True.
        do_loop_modelling (bool, optional): whether to do loop modelling. Defaults to False.
        charge_method (charge_methods, optional): charge method. Defaults to "bcc".
        lig_force_field (lig_force_fields, optional): ligand force field. Defaults to "gaff2".
        padding (float, optional): padding. Defaults to 1.0.
        system_name (str, optional): system name. Defaults to "complex".
        force_field (force_fields, optional): force field. Defaults to "ff14SB".

    """

    # input validation
    if padding < 0.5 or padding > 2:
        raise DeepOriginException("Padding must be greater than 0.5, and less than 2")

    database = _ensure_db_for_abfe()

    url = construct_resource_url(
        name=ABFE_DB,
        row_type="database",
    )

    print(f"Using row {row_id} in database at: {url}")

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
        "force": 2,
        "test_run": 0,
        "system": system_name,
        "include_ligands": include_ligands,
        "include_protein": include_protein,
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
            "columnId": output_column_name,
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
