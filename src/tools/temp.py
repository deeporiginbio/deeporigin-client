import json

import requests

url = "{{baseUrl}}/tools/{{orgFriendlyId}}"

payload = {
    "toolId": "deeporigin/md-suite-prep",
    "inputs": {
        "ligand": {
            "columnId": "_column:3clokepAhJikduTuFl6Wa",
            "rowId": "fma-1",
            "databaseId": "fep-mason-all",
        },
        "protein": {
            "columnId": "_column:A48UMj1CNibOJhe6GFaEL",
            "rowId": "fma-1",
            "databaseId": "fep-mason-all",
        },
        "force": 1,
        "test_run": 0,
        "system": "complex",
        "include_ligands": 1,
        "include_protein": 1,
        "sysprep_params": {
            "is_protein_protonated": True,
            "do_loop_modelling": False,
            "force_field": "ff14SB",
            "padding": 1,
            "keep_waters": True,
            "save_gmx_files": False,
            "is_lig_protonated": True,
            "charge_method": "bcc",
            "lig_force_field": "gaff2",
        },
    },
    "outputs": {
        "output_file": {
            "columnId": "_column:5rAPtk2oPl9psZZU3IeS0",
            "rowId": "fma-1",
            "databaseId": "fep-mason-all",
        }
    },
    "clusterId": {{cluster_id}},
}
headers = {
    "Content-Type": "application/json",
    "Authorization": "Bearer {{auth_token_edge}}",
}

response = requests.request("POST", url, headers=headers, data=payload)

print(response.text)
