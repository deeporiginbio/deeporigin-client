import json

import requests

url = "{{baseUrl}}/tools/{{orgFriendlyId}}"

payload = json.dumps(
    {
        "toolId": "md-suite-solvation",
        "inputs": {
            "input": {
                "columnId": "_column:5rAPtk2oPl9psZZU3IeS0",
                "rowId": "fma-2",
                "databaseId": "fep-mason-all",
            },
            "force": 1,
            "test_run": 0,
            "run_name": "annihilation",
            "system": "ligand_only",
            "softcore_alpha": 0.5,
            "annihilate": True,
            "em_solvent": True,
            "em_all": True,
            "nvt_heating_ns": 0.1,
            "npt_reduce_restraints_ns": 0.2,
            "repeats": 1,
            "steps": 300000,
            "threads": 0,
            "fep_windows": [
                {"coul_A": [1, 0.8, 0.6, 0.4, 0.2, 0]},
                {"vdw_A": [1, 0.9, 0.8, 0.7, 0.6, 0.5, 0.4, 0.3, 0.2, 0.1, 0]},
            ],
            "emeq_md_options": {
                "Δt": 0.004,
                "T": 298.15,
                "cutoff": 0.9,
                "fourier_spacing": 0.12,
                "hydrogen_mass": 2,
            },
            "prod_md_options": {
                "integrator": "BAOABIntegrator",
                "Δt": 0.004,
                "T": 298.15,
                "cutoff": 0.9,
                "fourier_spacing": 0.12,
                "hydrogen_mass": 2,
                "barostat": "MonteCarloBarostat",
                "barostat_exchange_interval": 500,
            },
        },
        "outputs": {
            "output_file": {
                "columnId": "_column:bjGXjeJQiCTXBBXxYjlO0",
                "rowId": "fma-2",
                "databaseId": "fep-mason-all",
            }
        },
        "clusterId": "d638d3cd-833f-476b-87fe-ca6d76e66049",
    }
)
headers = {
    "Content-Type": "application/json",
    "Authorization": "Bearer {{auth_token_edge}}",
}

response = requests.request("POST", url, headers=headers, data=payload)

print(response.text)
