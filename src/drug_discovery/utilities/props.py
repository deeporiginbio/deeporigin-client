# from typing import Optional, Dict, Any

# TODO: this needs to use the protonation and propertiy prediction functions from fast running tools

# def protonate(smiles: str, pH: float = 7.4, filter_percentage: float = 1):
#     """
#     Predicts the right protonation of a molecule at given pH value.

#     Parameters:
#     - entry: A single or multiple ligands represented as SMILES or Ligand instances.
#     - pH: pH value of the solvent for concentration calculation. Default is 7.4.
#     - filter_percentage: Percentage threshold for filtering low concentration states. Default is 1.

#     Returns:
#     - ProtonationReport: A ProtonationReport instance.
#     """
#     from deeporigin.src.properties import MolecularPropertiesClient

#     try:
#         client = MolecularPropertiesClient()
#         protonation_output = client.protonate(
#             entry=[smiles],
#             pH=pH,
#             filter_percentage=filter_percentage,
#             html_output=False,
#         )[0]
#         smiles = protonation_output["protonation"]["smiles_list"][0]

#         return smiles
#     except Exception as e:
#         return None


# def predict_properties(
#     smiles: str,
# ) -> Optional[Dict[str, Any]]:
#     """
#     Predicts specified properties for given ligands.

#     Parameters:
#     - smiles: A single SMILES strings.
#     - properties_to_predict: A dictionary specifying which properties to predict.
#                              Defaults to predicting hERG, logD, logP, logS, ames, and cyp.
#     - logger: An instance of Logger. If not provided, a default logger is initialized.

#     Returns:
#     - A dictionary containing the predicted properties if successful, else None.
#     """
#     properties_to_predict = {
#         "hERG": True,
#         "logD": True,
#         "logP": True,
#         "logS": True,
#         "ames": True,
#         "cyp": True,
#     }

#     request_body = {"smiles": smiles, "predict": properties_to_predict}

#     client = Client()
#     try:
#         response = client.post_request(endpoint="properties", data=request_body)
#     except Exception as e:
#         return None

#     if response.status_code == 200:
#         try:
#             response_data = response.json()
#             return response_data
#         except ValueError:
#             return None
#     else:
#         try:
#             error_message = response.json().get("msg", "Unknown error occurred.")
#         except ValueError:
#             pass
#         return None
