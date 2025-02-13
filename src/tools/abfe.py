"""this module contains various functions to run steps of an ABFE workflow"""

import os
from dataclasses import asdict, dataclass, field, fields
from pathlib import Path
from typing import List, Literal, Optional, Union

from beartype import beartype
from deeporigin import chemistry
from deeporigin.data_hub import api
from deeporigin.exceptions import DeepOriginException
from deeporigin.tools import run
from deeporigin.tools.toolkit import _ensure_columns, _ensure_database
from deeporigin.tools.utils import query_run_status
from deeporigin.utils.config import construct_resource_url
from deeporigin.utils.notebook import render_mermaid
from IPython.display import HTML, display

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


@beartype
def _ensure_db_for_abfe() -> dict:
    """ensure that there is a database for FEP on Data Hub"""

    database = _ensure_database(ABFE_DB)

    try:
        api.add_smiles_column(
            name="Ligand",
            database_id=ABFE_DB,
        )
    except Exception:
        pass

    required_columns = [
        dict(name="FEP ΔG (kcal/mol)", type="float"),
        dict(name="ligand_file", type="file"),
        dict(name="protein_file", type="file"),
        dict(name="complex_prep_output", type="file"),
        dict(name="ligand_prep_output", type="file"),
        dict(name="emeq_output", type="file"),
        dict(name="solvation_output", type="file"),
        dict(name="md_output", type="file"),
        dict(name="abfe_output", type="file"),
        dict(name="end_to_end_output", type="file"),
    ]

    database = _ensure_columns(
        database=database,
        required_columns=required_columns,
    )

    return database


@dataclass
class MDParams:
    integrator: integrators = "BAOABIntegrator"
    Δt: float = 0.004
    T: float = 298.15
    cutoff: float = 0.9
    fourier_spacing: float = 0.12
    hydrogen_mass: int = 2
    barostat: str = "MonteCarloBarostat"
    barostat_exchange_interval: int = 500


@dataclass
class SystemPrepParams:
    charge_method: charge_methods = "bcc"
    do_loop_modelling: bool = False
    force_field: force_fields = "ff14SB"
    is_lig_protonated: bool = True
    is_protein_protonated: bool = True
    keep_waters: bool = False
    lig_force_field: ligand_force_fields = "gaff2"
    padding: float = 1.0  # nm
    save_gmx_files: bool = False


@dataclass
class Ligand:
    file: Union[str, Path]
    smiles_string: Optional[str] = None
    is_protonated: Optional[bool] = False

    def __post_init__(self):
        if self.smiles_string is None:
            self.smiles_string = chemistry.sdf_to_smiles(self.file)

        if self.is_protonated is None:
            self.is_protonated = chemistry.is_ligand_protonated(self.file)

    def _repr_pretty_(self, p, cycle):
        if cycle:
            p.text("Ligand(...)")
        else:
            p.text("Ligand(")

            with p.group(2, "\n  ", "\n"):
                all_fields = fields(self)
                for idx, field in enumerate(all_fields):
                    value = getattr(self, field.name)
                    p.text(f"{field.name}: {value!r}")
                    # Only add a breakable if this isn't the last field.
                    if idx < len(all_fields) - 1:
                        p.breakable()
            p.text(")")


@dataclass
class Protein:
    file: str | Path


@dataclass
class ABFE:
    """ABFE class that can work with one protein and many ligands"""

    ligands: List[Ligand]
    protein: Protein

    row_ids: Optional[list[str]] = None

    _job_ids: dict = field(default_factory=dict)
    _status: dict = field(default_factory=dict)

    @classmethod
    def from_dir(cls, directory: str) -> "ABFE":
        """initialize an ABFE class given some files in a directory"""

        sdf_files = sorted(
            [
                os.path.join(directory, f)
                for f in os.listdir(directory)
                if f.lower().endswith(".sdf")
            ]
        )
        ligands = [Ligand(sdf_file) for sdf_file in sdf_files]

        pdb_files = [
            os.path.join(directory, f)
            for f in os.listdir(directory)
            if f.lower().endswith(".pdb")
        ]

        if len(pdb_files) != 1:
            raise ValueError(
                f"Expected exactly one PDB file in the directory, but found {len(pdb_files)}."
            )
        protein_file = pdb_files[0]
        protein = Protein(protein_file)

        return cls(ligands=ligands, protein=protein)

    def init(self):
        """Initialize an ABFE run. Upload ligand(s) and protein files to Data Hub."""

        database = _ensure_db_for_abfe()

        smiles_column_id = [col.id for col in database.cols if col.name == "Ligand"][0]

        url = construct_resource_url(
            name=ABFE_DB,
            row_type="database",
        )

        print(f"Using database at: {url}")
        print("🧬 Uploading files to database...")

        self.row_ids = []

        # now upload the protein
        protein_file = api.upload_file(file_path=self.protein.file)

        for ligand in self.ligands:
            # for each ligand, upload the ligand to a new row
            response = api.upload_file_to_new_database_row(
                database_id=ABFE_DB,
                column_id="ligand_file",
                file_path=ligand.file,
            )

            row_id = response.rows[0].hid

            self.row_ids.append(row_id)

            api.set_cell_data(
                ligand.smiles_string,
                database_id=ABFE_DB,
                row_id=row_id,
                column_id=smiles_column_id,
            )

            api.assign_files_to_cell(
                file_ids=[protein_file.id],
                database_id=ABFE_DB,
                column_id="protein_file",
                row_id=row_id,
            )

        print(f"🧬 Files uploaded to row {row_id}.")

        # self._status["init"] = "Succeeded"

    def complex_prep(
        self,
        params: SystemPrepParams = SystemPrepParams(),
    ):
        """Run complex prep on a ligand and protein pair, that exist as files on a row in the ABFE database. For this to work, the complex prep step must have been run first."""

        self._job_ids["complex_prep"] = _prep(
            name="complex",
            row_id=self.row_id,
            output_column_name="complex_prep_output",
            sysprep_params=params,
            include_protein=1,
            include_ligands=1,
        )
        self._status["complex_prep"] = "NotStarted"

    @beartype
    def ligand_prep(
        self,
        params: SystemPrepParams = SystemPrepParams(),
    ):
        self._job_ids["ligand_prep"] = _prep(
            name="complex",
            row_id=self.row_id,
            output_column_name="ligand_prep_output",
            sysprep_params=params,
            include_protein=0,
            include_ligands=1,
        )
        self._status["ligand_prep"] = "NotStarted"

    @beartype
    def emeq(
        self,
        *,
        system_name: str = "complex",
        nvt_heating_ns: float = 0.1,
        npt_reduce_restraints_ns: float = 0.2,
        params: MDParams = MDParams(),
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

        params = asdict(params)

        tool_key = "deeporigin.md-suite-emeq"

        database = _ensure_database(ABFE_DB)

        inputs = {
            "input": {
                "columnId": "complex_prep_output",
                "rowId": self.row_id,
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
            "emeq_md_options": params,
        }

        outputs = {
            "output_file": {
                "columnId": "emeq_output",
                "rowId": self.row_id,
                "databaseId": database.hid,
            }
        }

        self._job_ids["emeq"] = run._process_job(
            inputs=inputs,
            outputs=outputs,
            tool_key=tool_key,
            cols=database.cols,
        )

        self._status["emeq"] = "NotStarted"

    def status(self):
        node_statuses = {
            "init": "NotStarted",
            "complex_prep": "NotStarted",
            "ligand_prep": "NotStarted",
            "solvation_FEP": "NotStarted",
            "simple_MD": "NotStarted",
            "binding_FEP": "NotStarted",
            "DeltaG": "NotStarted",
        }

        keys = self._job_ids.keys()

        for key in keys:
            if self._status[key] == "Succeeded" or self._status[key] == "Failed":
                continue

            # IN NON-terminal state. update
            self._status[key] = query_run_status(self._job_ids[key])

        node_statuses.update(self._status)

        render_mermaid_with_statuses(node_statuses)


@beartype
def _prep(
    *,
    name: str,
    row_id: str,
    output_column_name: str,
    sysprep_params: SystemPrepParams,
    include_protein: int = 1,
    include_ligands: int = 1,
) -> str:
    """Function to prepare uploaded Ligand and protein files using Deep Origin MDSuite. Use this function to run system prep on a ligand and protein pair, that exist as files on a row in the ABFE database."""

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
        "system": name,
        "include_ligands": include_ligands,
        "include_protein": include_protein,
        "sysprep_params": asdict(sysprep_params),
    }

    outputs = {
        "output_file": {
            "columnId": output_column_name,
            "rowId": row_id,
            "databaseId": database.hid,
        }
    }

    job_id = run._process_job(
        inputs=inputs,
        outputs=outputs,
        tool_key=tool_key,
        cols=database.cols,
    )

    return job_id


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


@beartype
def render_mermaid_with_statuses(
    node_status: dict[str, str] = None,
) -> None:
    """
    Render a Mermaid diagram where each node is drawn as a rounded rectangle
    with a color indicating its status.

    Any node not specified in the node_status dict will default to "notStarted".

    """
    # Define the fixed nodes in the diagram.
    nodes = [
        "init",
        "emeq",
        "complex_prep",
        "ligand_prep",
        "solvation_FEP",
        "simple_MD",
        "binding_FEP",
        "DeltaG",
    ]

    # Use an empty dictionary if none is provided.
    if node_status is None:
        node_status = {}

    # Build node definitions. For each node, use the provided status or default to "notStarted".
    node_defs = ""
    for node in nodes:
        status = node_status.get(node, "NotStarted")
        node_defs += f"    {node}({node}):::{status};\n"

    # Define the fixed edges of the diagram.
    edges = """
    init --> complex_prep;
    init --> ligand_prep;
    ligand_prep ----> solvation_FEP;
    solvation_FEP --> DeltaG;
    complex_prep --> emeq --> simple_MD --> binding_FEP --> DeltaG;
    """

    # Build the complete Mermaid diagram definition.
    mermaid_code = f"""
graph LR;
    %% Define styles for statuses:
    classDef NotStarted fill:#cccccc,stroke:#333,stroke-width:2px;
    classDef Queued fill:#cccccc,stroke:#222,stroke-width:2px;
    classDef Succeeded    fill:#90ee90,stroke:#333,stroke-width:2px;
    classDef Running       fill:#87CEFA,stroke:#333,stroke-width:2px;
    classDef failed     fill:#ff7f7f,stroke:#333,stroke-width:2px;

{node_defs}
{edges}
    """

    # Render the diagram using your helper function.
    render_mermaid(mermaid_code)

    # Define HTML for the legend. Each status is displayed as a colored span.
    legend_html = """
    <div style="margin-top: 20px; font-family: sans-serif;">
      <span style="background-color:#cccccc; color: black; padding:2px 4px; margin: 0 8px;">NotStarted</span>
      <span style="background-color:#90ee90; color: black; padding:2px 4px; margin: 0 8px;">Suceedeed</span>
      <span style="background-color:#87CEFA; color: black; padding:2px 4px; margin: 0 8px;">Running</span>
      <span style="background-color:#ff7f7f; color: black; padding:2px 4px; margin: 0 8px;">Failed</span>
    </div>
    """
    # Display the legend below the Mermaid diagram.
    display(HTML(legend_html))
