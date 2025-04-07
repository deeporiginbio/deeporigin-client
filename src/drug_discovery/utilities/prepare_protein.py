# from typing import List

# from deeporigin.src.client import Client
# from deeporigin.src.utilities.logging import DEFAULT_LOGGER


# def prepare(
#     protein_path: str,
#     protein_pdb_id: str = "",
#     protein_extension: str = "pdb",
#     metal_resnames: List[str] = None,
#     cofactor_resnames: List[str] = None,
#     model_loops: bool = False,
# ):
#     client = Client()
#     with open(protein_path, "r") as e:
#         data = e.read()

#     metal_resnames = ",".join(metal_resnames) if metal_resnames else ""
#     cofactor_resnames = ",".join(cofactor_resnames) if cofactor_resnames else ""

#     payload = {
#         "content": data,
#         "pdb_id": protein_pdb_id,
#         "extension": protein_extension,
#         "metals": metal_resnames,
#         "cofactors": cofactor_resnames,
#         "model_loops": model_loops,
#     }

#     response = client.post_request(
#         logger=DEFAULT_LOGGER,
#         endpoint="docking/prepare",
#         data=payload,
#     )
#     pdb_content = ""
#     if response.status_code == 200:
#         data = response.json()
#         DEFAULT_LOGGER.log_info(data["msg"])

#         pdb_content = data["pdb_content"]
#     else:
#         msg = data.get("msg", "")
#         if not msg:
#             msg = "Failed to prepare"
#         DEFAULT_LOGGER.log_error(msg)

#     return {
#         "prepared_protein_content": pdb_content,
#         "raw_protein_path": protein_path,
#         "protein_extension": protein_extension,
#         "protein_pdb_id": protein_pdb_id,
#     }
