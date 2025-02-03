"""this module contains various functions to run steps of an ABFE workflow"""

from typing import Literal

from beartype import beartype
from deeporigin.data_hub import api
from deeporigin.exceptions import DeepOriginException
from deeporigin.tools import run
from deeporigin.tools.toolkit import _ensure_columns, _ensure_database
from deeporigin.utils.config import construct_resource_url

# constants and types
ABFE_DB = "ABFE"

charge_methods = Literal["gas", "bcc"]
"""Available charge methods"""

ligand_force_fields = Literal["gaff", "gaff2", "openff"]
"""Available ligand force fields. `gaff` is General Amber Force Field, gaff2 is an updated version of gaff. """

force_fields = Literal["ff14SB", "ff99SB-ildn"]
"""Available force fields"""

integrators = Literal[
    "BAOABIntegrator",
    "LangevinIntegrator",
    "SteepestDescentMinimizer",
    "LeapFrogIntegrator",
    "SoluteSolventSplittingIntegrator",
]
"""Integrator available for simulation"""


prod_md_defaults = {
    "integrator": "BAOABIntegrator",
    "Δt": 0.004,
    "T": 298.15,
    "cutoff": 0.9,
    "fourier_spacing": 0.12,
    "hydrogen_mass": 2,
    "barostat": "MonteCarloBarostat",
    "barostat_exchange_interval": 500,
}


emeq_md_defaults = {
    "Δt": 0.004,
    "T": 298.15,
    "cutoff": 0.9,
    "fourier_spacing": 0.12,
    "hydrogen_mass": 2,
}


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
        dict(name="solvation_output", type="file"),
        dict(name="md_output", type="file"),
        dict(name="abfe_output", type="file"),
    ]

    database = _ensure_database(ABFE_DB)

    database = _ensure_columns(
        database=database,
        required_columns=required_columns,
    )

    return database


@beartype
def emeq(
    row_id: str,
    *,
    system_name: str = "complex",
    fourier_spacing: float = 0.12,
    hydrogen_mass: int = 2,
    cutoff: float = 0.9,
    T: float = 298.15,
    Δt: float = 0.004,
    npt_reduce_restraints_ns: float = 0.2,
    nvt_heating_ns: float = 0.1,
) -> None:
    """Run emeq on a ligand and protein pair, that exist as files on a row in the ABFE database. For this to work, the complex prep step must have been run first.

    Args:
        row_id (str): row id that contains the ligand and protein files.
        system_name (str, optional): name of the system. Defaults to "complex". This name can be anything.
        fourier_spacing (float, optional): spacing of the fourier grid. Defaults to 0.12.
        hydrogen_mass (int, optional): hydrogen mass. Defaults to 2.
        cutoff (float, optional): cutoff. Defaults to 0.9.
        T (float, optional): temperature. Defaults to 298.15.
        Δt (float, optional): time step. Defaults to 0.004.
        npt_reduce_restraints_ns (float, optional): time step. Defaults to 0.2.
        nvt_heating_ns (float, optional): time step. Defaults to 0.1.


    """

    kwargs = locals()

    tool_key = "deeporigin.md-suite-emeq"

    database = _ensure_database(ABFE_DB)

    emeq_md_options = {
        k: kwargs[k] if k in kwargs else v for k, v in emeq_md_defaults.items()
    }

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
        "nvt_heating_ns": nvt_heating_ns,
        "npt_reduce_restraints_ns": npt_reduce_restraints_ns,
        "from_run": "__USE_SYSTEM",
        "emeq_md_options": emeq_md_options,
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
        tool_key=tool_key,
        cols=database.cols,
    )


@beartype
def complex_prep(
    row_id: str,
    *,
    keep_waters: bool = True,
    save_gmx_files: bool = False,
    is_lig_protonated: bool = True,
    is_protein_protonated: bool = True,
    do_loop_modelling: bool = False,
    charge_method: charge_methods = "bcc",
    lig_force_field: ligand_force_fields = "gaff2",
    padding: float = 1.0,
    system_name: str = "complex",
    force_field: force_fields = "ff14SB",
) -> None:
    """Function to prepare uploaded Ligand and protein files using Deep Origin MDSuite. Use this function to run system prep on a ligand and protein pair, that exist as files on a row in the ABFE database.

    Args:
        row_id (str): row id that contains the ligand and protein files.
        keep_waters (bool, optional): whether to keep water molecules in the system. Defaults to True.
        save_gmx_files (bool, optional): whether to save gmx files. Defaults to False.
        is_lig_protonated (bool, optional): whether the ligand is protonated. Defaults to True.
        is_protein_protonated (bool, optional): whether the protein is protonated. Defaults to True.
        do_loop_modelling (bool, optional): whether to do loop modelling. Defaults to False.
        charge_method (str, optional): method to use for charge assignment. Defaults to "bcc".
        lig_force_field (str, optional): ligand force field. Defaults to "gaff2".
        padding (float, optional): padding to use. Defaults to 1.0.
        force_field (str, optional): force field. Defaults to "ff14SB".
        system_name (str, optional): name of the system. Defaults to "complex". This name can be anything.


    """

    kwargs = locals()

    kwargs["include_protein"] = 1
    kwargs["include_ligands"] = 1
    kwargs["output_column_name"] = "complex_prep_output"

    _prep(**kwargs)


@beartype
def solvation_prep(
    row_id: str,
    *,
    is_lig_protonated: bool = True,
    charge_method: charge_methods = "bcc",
    lig_force_field: ligand_force_fields = "gaff2",
    system_name: str = "ligand_prep",
    force_field: force_fields = "ff14SB",
) -> None:
    """Function to prepare uploaded Ligand and protein files using Deep Origin MDSuite. Use this function to run system prep on a ligand and protein pair, that exist as files on a row in the ABFE database.

    Args:
        row_id (str): row id that contains the ligand and protein files.
        is_lig_protonated (bool, optional): whether the ligand is protonated. Defaults to True.
        charge_method (str, optional): method to use for charge assignment. Defaults to "bcc".
        lig_force_field (str, optional): ligand force field. Defaults to "gaff2".
        system_name (str, optional): name of the system. Defaults to "ligand_prep". This name can be anything.
        force_field (str, optional): force field. Defaults to "ff14SB".


    """
    kwargs = locals()

    kwargs["include_protein"] = 0
    kwargs["include_ligands"] = 1
    kwargs["output_column_name"] = "solvation_prep_output"

    _prep(**kwargs)


@beartype
def _prep(
    row_id: str,
    *,
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
    lig_force_field: ligand_force_fields = "gaff2",
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
        lig_force_field (ligand_force_fields, optional): ligand force field. Defaults to "gaff2".
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

    tool_key = "deeporigin.md-suite-prep"

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
        tool_key=tool_key,
        cols=database.cols,
    )


def solvation_fep(
    row_id: str,
    *,
    integrator: integrators = "BAOABIntegrator",
    softcore_alpha: float = 0.5,
    steps: int = 300000,
    repeats: int = 1,
    threads: int = 0,
    annihilate: bool = True,
    hydrogen_mass: int = 2,
    T: float = 298.15,
) -> None:
    """Run a solvation simulation


    Args:
        row_id (str): row id of the ligand and protein files.
        integrator (integrators, optional): integrator. Defaults to "BAOABIntegrator".
        softcore_alpha (float, optional): softcore alpha. Defaults to 0.5.
        steps (int, optional): The number of steps to run the simulation for the prod step.
        repeats (int, optional): The number of repeats for prod step.
        threads (int, optional): The number of threads per worker. By default the number of threads will be determined by the number of windows and available cores on the CPU. Defaults to 0.
        annihilate (bool, optional): Whether to annihilate the ligand or decouple it. Defaults to True.
    """

    kwargs = locals()

    if softcore_alpha < 0.0 or softcore_alpha > 1.0:
        raise DeepOriginException("softcore_alpha must be between 0.0 and 1.0")

    database = _ensure_db_for_abfe()

    url = construct_resource_url(
        name=ABFE_DB,
        row_type="database",
    )

    print(f"Using row {row_id} in database at: {url}")

    tool_key = "deeporigin.md-suite-solvation"

    prod_md_options = {
        k: kwargs[k] if k in kwargs else v for k, v in prod_md_defaults.items()
    }
    emeq_md_options = {
        k: kwargs[k] if k in kwargs else v for k, v in emeq_md_defaults.items()
    }

    inputs = {
        "input": {
            "columnId": "solvation_prep_output",
            "rowId": row_id,
            "databaseId": database.hid,
        },
        "force": 1,
        "test_run": 0,
        "run_name": "annihilation",
        "system": "ligand_only",
        "softcore_alpha": softcore_alpha,
        "annihilate": annihilate,
        "em_solvent": True,
        "em_all": True,
        "nvt_heating_ns": 0.1,
        "npt_reduce_restraints_ns": 0.2,
        "repeats": repeats,
        "steps": steps,
        "threads": threads,
        "fep_windows": [
            {"coul_A": [1, 0.8, 0.6, 0.4, 0.2, 0]},
            {"vdw_A": [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0]},
        ],
        "emeq_md_options": emeq_md_options,
        "prod_md_options": prod_md_options,
    }

    outputs = {
        "output_file": {
            "columnId": "solvation_output",
            "rowId": row_id,
            "databaseId": database.hid,
        }
    }

    run._process_job(
        inputs=inputs,
        outputs=outputs,
        tool_key=tool_key,
        cols=database.cols,
    )


def simple_md(
    row_id: str,
    *,
    integrator: integrators = "BAOABIntegrator",
    run_name: str = "complex_1ns_md",
    steps: int = 250000,
    threads: int = 0,
    Δt: float = 0.004,
    temperature: float = 298.15,
    cutoff: float = 0.9,
    fourier_spacing: float = 0.12,
):
    """Run a simple MD simulation

    Args:
        row_id (str): row id of the ligand and protein files.
        integrator (integrators, optional): integrator. Defaults to "BAOABIntegrator".
        run_name (str, optional): run name. Defaults to "complex_1ns_md".
        steps (int, optional): The number of steps to run the simulation for the prod step.
        threads (int, optional): The number of threads per worker. By default the number of threads will be determined by the number of windows and available cores on the CPU. Defaults to 0.
        Δt (float, optional): The time step in femtoseconds. Defaults to 0.004.
        temperature (float, optional): The temperature in kelvin. Defaults to 298.15.
        cutoff (float, optional): The cutoff distance in angstroms. Defaults to 0.9.
        fourier_spacing (float, optional): The Fourier spacing in femtoseconds. Defaults to 0.12.
    """
    kwargs = locals()
    database = _ensure_db_for_abfe()

    url = construct_resource_url(
        name=ABFE_DB,
        row_type="database",
    )

    md_options = {
        k: kwargs[k] if k in kwargs else v for k, v in prod_md_defaults.items()
    }

    print(f"Using row {row_id} in database at: {url}")

    tool_key = "deeporigin.md-suite-md"
    inputs = {
        "input": {
            "columnId": "emeq_output",
            "rowId": row_id,
            "databaseId": ABFE_DB,
        },
        "force": 1,
        "test_run": 0,
        "run_name": run_name,
        "from_run": "test-run",
        "steps": steps,
        "threads": threads,
        "system": "__USE_FROM_RUN",
        "md_options": md_options,
    }

    outputs = {
        "output_file": {
            "columnId": "md_output",
            "rowId": row_id,
            "databaseId": database.hid,
        }
    }

    run._process_job(
        inputs=inputs,
        outputs=outputs,
        tool_key=tool_key,
        cols=database.cols,
    )


def binding_fep(
    row_id,
    *,
    run_name: str = "annihilation_fep",
    softcore_alpha: float = 0.5,
    annihilate: bool = True,
    em_solvent: bool = True,
    em_all: bool = True,
    nvt_heating_ns: int = 1,
    npt_reduce_restraints_ns: int = 2,
    steps: int = 1250000,
    repeats: int = 1,
    threads: int = 0,
    integrator: integrators = "BAOABIntegrator",
    Δt: float = 0.004,
    T: float = 298.15,
):
    """Run an ABFE simulation

    Args:
        row_id (str): row id of the ligand and protein files.
        run_name (str, optional): run name. Defaults to "annihilation_fep".
        softcore_alpha (float, optional): softcore alpha. Defaults to 0.5.
        annihilate (bool, optional): Whether to annihilate the ligand or decouple it. Defaults to True.
        em_solvent (bool, optional): Whether to use em solvent. Defaults to True.
        em_all (bool, optional): Whether to use em all. Defaults to True.
        nvt_heating_ns (int, optional): The number of nvt heating steps. Defaults to 1.
        npt_reduce_restraints_ns (int, optional): The number of npt reduce restraints steps. Defaults to 2.
        steps (int, optional): The number of steps to run the simulation for the prod step.
        threads (int, optional): The number of threads per worker. By default the number of threads will be determined by the number of windows and available cores on the CPU. Defaults to 0.
        integrator (integrators, optional): integrator. Defaults to "BAOABIntegrator".
        Δt (float, optional): The time step in femtoseconds. Defaults to 0.004.
        T (float, optional): The temperature in kelvin. Defaults to 298.15.

    """

    kwargs = locals()
    database = _ensure_db_for_abfe()

    url = construct_resource_url(
        name=ABFE_DB,
        row_type="database",
    )

    print(f"Using row {row_id} in database at: {url}")

    emeq_md_options = {
        k: kwargs[k] if k in kwargs else v for k, v in emeq_md_defaults.items()
    }
    prod_md_options = {
        k: kwargs[k] if k in kwargs else v for k, v in prod_md_defaults.items()
    }

    tool_key = "deeporigin.md-suite-abfe"
    inputs = {
        "input": {
            "columnId": "md_output",
            "rowId": row_id,
            "databaseId": ABFE_DB,
        },
        "force": 1,
        "test_run": 0,
        "run_name": run_name,
        "boresch_run_name": "complex_1ns_md",
        "system": "complex",
        "softcore_alpha": softcore_alpha,
        "annihilate": annihilate,
        "em_solvent": em_solvent,
        "em_all": em_all,
        "nvt_heating_ns": nvt_heating_ns,
        "npt_reduce_restraints_ns": npt_reduce_restraints_ns,
        "repeats": repeats,
        "steps": steps,
        "threads": threads,
        "fep_windows": [
            {"restraints_A": [0, 0.01, 0.025, 0.05, 0.1, 0.35, 0.5, 0.75, 1]},
            {"coul_A": [1, 0.8, 0.6, 0.4, 0.2, 0]},
            {"vdw_A": [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0]},
        ],
        "emeq_md_options": emeq_md_options,
        "prod_md_options": prod_md_options,
    }

    outputs = {
        "output_file": {
            "columnId": "abfe_output",
            "rowId": row_id,
            "databaseId": database.hid,
        }
    }

    run._process_job(
        inputs=inputs,
        outputs=outputs,
        tool_key=tool_key,
        cols=database.cols,
    )
