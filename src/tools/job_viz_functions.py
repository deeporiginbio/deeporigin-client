"""helper module to parse progress and render progress for first party tools"""

import json
import os

from beartype import beartype


def _abfe_parse_progress(job) -> dict:
    """parse progress from a ABFE job"""

    steps = [
        "init",
        "complex",
        "ligand",
        "simple_md",
        "solvation",
        "binding",
        "delta_g",
    ]

    if job._progress_report is None:
        return dict.fromkeys(steps, "NotStarted")

    try:
        data = job._progress_report

        if data is None:
            progress = dict.fromkeys(steps, "NotStarted")
            progress["init"] = "Running"
            return progress
        else:
            data = json.loads(data)

        if "cmd" in data and data["cmd"] == "FEP Results":
            return dict.fromkeys(steps, "Succeeded")

        if "status" in data and data["status"] == "Initiating":
            progress = dict.fromkeys(steps, "NotStarted")
            progress["init"] = "Running"
            return progress

        status_value = job._status

        # If the overall status is Succeeded, return a dictionary with every key set to "Succeeded".
        if status_value == "Succeeded":
            return dict.fromkeys(steps, "Succeeded")

        current_step = data["run_name"]

        # Validate the input step
        if current_step not in steps:
            raise ValueError(
                f"Invalid process step provided: {current_step}. "
                f"Valid steps are: {', '.join(steps)}."
            )

        progress = {}
        for step in steps:
            if step == current_step:
                progress[step] = "Running"
                # Once we hit the current step, stop processing further steps.
                break
            else:
                progress[step] = "Succeeded"

        # if the job failed, mark the step that is running as failed
        if job._status == "Failed":
            progress[current_step] = "Failed"

    except Exception:
        progress = dict.fromkeys(steps, "Indeterminate")
        progress["init"] = "Indeterminate"

    return progress


@beartype
def _viz_func_rbfe(job) -> str:
    """
    Render HTML for a Mermaid diagram where each node is drawn as a rounded rectangle
    with a color indicating its status.
    """
    import json

    # For single job, we have single metadata and progress report
    metadata = job._metadata
    report = job._progress_report

    ligand1 = (
        metadata.get("ligand1_file", "Unknown ligand") if metadata else "Unknown ligand"
    )
    ligand2 = (
        metadata.get("ligand2_file", "Unknown ligand") if metadata else "Unknown ligand"
    )

    if report is None:
        step = ""
        sub_step = ""
    else:
        data = json.loads(report)
        step = data.get("cmd", "")
        sub_step = data.get("sub_step", "")

    import pandas as pd

    df = pd.DataFrame(
        {
            "ligand1": [ligand1],
            "ligand2": [ligand2],
            "steps": [step],
            "sub_steps": [sub_step],
        }
    )
    return df.to_html()


@beartype
def _viz_func_abfe(job) -> str:
    """
    Render HTML for a Mermaid diagram where each node is drawn as a rounded rectangle
    with a color indicating its status.

    Any node not specified in the node_status dict will default to "notStarted".
    """

    from deeporigin.utils.notebook import mermaid_to_html

    statuses = _abfe_parse_progress(job)

    # Define the fixed nodes in the diagram
    nodes = [
        "init(Init)",
        "complex(Complex Prep)",
        "ligand(Ligand Prep)",
        "solvation(Solvation FEP)",
        "simple_md(Simple MD)",
        "binding(Binding FEP)",
        "delta_g(ΔG)",
    ]

    # Build node definitions. For each node, use the provided status or default to "notStarted".
    node_defs = ""
    for node in nodes:
        label = node.split("(")[0]
        status = statuses.get(label, "NotStarted")
        node_defs += f"    {node}:::{status};\n"

    # Define the fixed edges of the diagram.
    edges = """
    init --> complex;
    init --> ligand;
    ligand ----> solvation;
    solvation --> delta_g;
    complex ---> simple_md --> binding -->delta_g;
    """

    # Build the complete Mermaid diagram definition.
    mermaid_code = f"""
graph LR;
    %% Define styles for statuses:
    classDef NotStarted   fill:#cccccc,stroke:#333,stroke-width:2px;
    classDef Queued    fill:#cccccc,stroke:#222,stroke-width:2px;
    classDef Succeeded   fill:#90ee90,stroke:#333,stroke-width:2px;
    classDef Running      fill:#87CEFA,stroke:#333,stroke-width:2px;
    classDef Failed    fill:#ff7f7f,stroke:#333,stroke-width:2px;

{node_defs}
{edges}
    """

    # Render the diagram using your helper function.
    mermaid_html = mermaid_to_html(mermaid_code)

    # Define HTML for the legend. Each status is displayed asa colored span.
    legend_html = """
    <div style="margin-top: 20px; font-family: sans-serif;">
        <span style="background-color:#cccccc; color: black;padding:2px 4px; margin: 0 8px;">NotStarted</span>
        <span style="background-color:#90ee90; color: black;padding:2px 4px; margin: 0 8px;">Succeeded</span>
        <span style="background-color:#87CEFA; color: black;padding:2px 4px; margin: 0 8px;">Running</span>
        <span style="background-color:#ff7f7f; color: black;padding:2px 4px; margin: 0 8px;">Failed</span>
    </div>
    """
    # Display the legend below the Mermaid diagram.
    return mermaid_html + legend_html


def _viz_func_docking(job) -> str:
    """Render progress visualization for a docking job."""

    data = job._progress_report

    total_ligands = len(job._inputs["smiles_list"]) if job._inputs else 0
    total_docked = 0
    total_failed = 0

    if data is not None:
        total_docked += data.count("ligand docked")
        total_failed += data.count("ligand failed")

    total_running_time = job._get_running_time() or 0
    speed = total_docked / total_running_time if total_running_time > 0 else 0

    from deeporigin.utils.notebook import render_progress_bar

    return render_progress_bar(
        completed=total_docked,
        total=total_ligands,
        failed=total_failed,
        title="Docking Progress",
        body_text=f"Average speed: {speed:.2f} dockings/minute",
    )


@beartype
def _name_func_docking(job) -> str:
    """Generate a name for a docking job."""

    unique_smiles = set()
    if job._inputs:
        unique_smiles.update(job._inputs["smiles_list"])
    num_ligands = len(unique_smiles)

    protein_file = (
        os.path.basename(job._metadata["protein_file"])
        if job._metadata
        else "Unknown protein"
    )

    return f"Docking <code>{protein_file}</code> to {num_ligands} ligands."


@beartype
def _name_func_abfe(job) -> str:
    """utility function to name a job using inputs to that job"""
    try:
        return f"ABFE run using <code>{job._metadata['protein_name']}</code> and <code>{job._metadata['ligand_name']}</code>"
    except Exception:
        return "ABFE run"


@beartype
def _name_func_rbfe(job) -> str:
    """utility function to name a job using inputs to that job"""

    try:
        # For single job, we always have single ligand pair
        return f"RBFE run using <code>{job._metadata['protein_file']}</code> and <code>{job._metadata['ligand_file']}</code>"
    except Exception:
        return "RBFE run"
