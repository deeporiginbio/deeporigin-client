# ┌▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀▀┐
# █ Results and Reports                                                                                                      █
# █                                                                                                                          █
# █ Description: This module provides comprehensive tools for molecular docking analysis and visualization. It               █
# █ includes classes for handling docking results, generating detailed reports, calculating molecular properties,            █
# █ and visualizing protein-ligand interactions. The primary classes and their functionalities are as follows:               █
# █                                                                                                                          █
# █ - **DockingResult**: Represents the outcome of a docking simulation, managing ligands, calculating RMSDs,                █
# █   exporting results to SDF files, and generating dataframes for analysis.                                                █
# █                                                                                                                          █
# █ - **PocketFinderReport**: Handles the identification and ensemble creation of binding pockets within proteins,           █
# █   and provides visualization capabilities for the detected pockets.                                                      █
# █                                                                                                                          █
# █ - **MolPropsReport**: Aggregates and formats molecular properties predictions, facilitating easy data analysis.          █
# █                                                                                                                          █
# █ - **ProtonationReport**: Manages protonation state predictions of molecules at specified pH levels, including            █
# █   visualization of concentration curves and molecular structures.                                                        █
# █                                                                                                                          █
# █ - **PainsReport**: Detects and highlights PAINS (Pan Assay Interference Compounds) substructures within molecules,       █
# █   providing visual representations and detailed reports.                                                                 █
# █                                                                                                                          █
# █ - **DockingReport**: Consolidates multiple DockingResult instances, enabling comprehensive reporting, data               █
# █   exportation, and visualization of aggregated docking data.                                                             █
# █                                                                                                                          █
# █ The module leverages various libraries such as RDKit for cheminformatics, Plotly for interactive visualizations,         █
# █ and Pandas for data manipulation. It also integrates with Jupyter for seamless visualization within notebooks.           █
# └▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄▄┘

import base64
import copy
import os
import subprocess
import tempfile
from copy import deepcopy
from datetime import datetime
from io import BytesIO
from pathlib import Path
from typing import List, Optional

import numpy as np
import pandas as pd
import plotly.graph_objs as go
import plotly.subplots as sp
import prolif as plf
from deeporigin_molstar import DockingViewer
from IPython.display import HTML, display
from prolif.interactions import VdWContact
from rdkit import Chem
from rdkit.Chem import AllChem, Draw

from ..utilities.alignments import save_bounding_box
from ..utilities.logging import DEFAULT_LOGGER
from ..utilities.utils import move_file_with_extension, remove_file
from ..utilities.visualize import jupyter_visualization
from .ligand import Ligand
from .protein import Protein


class DockingResult:
    def __init__(
        self,
        protein: Protein,
        smiles: Optional[str] = None,
        file_path: Optional[str] = None,
        successful: Optional[bool] = True,
    ):
        self.protein = protein
        self.ligands: List[Ligand] = []
        self.rmsds: Optional[List[float]] = None
        self.top_ligand: Optional[Ligand] = None
        self.smiles = smiles
        self.successful = successful
        self.file_path = file_path

    def add_ligand(self, ligand: Ligand):
        self.ligands.append(ligand)

    def calculate_rmsds_from_crystal(self, crystal_ligand: Ligand | str):
        if isinstance(crystal_ligand, str):
            crystal_ligand = Ligand(file_path=str(crystal_ligand))

        try:
            result = subprocess.run(
                ["obrms", self.file_path, crystal_ligand.file_path],
                capture_output=True,
                text=True,
            )
            self.rmsds = [
                float(rmsd.split()[-1]) for rmsd in result.stdout.split("\n") if rmsd
            ]
            return self.rmsds
        except Exception as e:
            raise SystemError(f"Failed to calculate RMSD values: {e}")

    def _to_sdf(self, safe=True, sdf_file_path=None):
        if not self.ligands:
            return None

        if not sdf_file_path:
            sdf_file_path = os.path.join(
                os.path.dirname(self.file_path or tempfile.gettempdir()),
                f"{self.protein.name}_docking_result_{self.smiles}.sdf",
            )

        if safe and os.path.isfile(sdf_file_path):
            move_file_with_extension(sdf_file_path, "sdf")
        else:
            remove_file(sdf_file_path)

        writer = Chem.SDWriter(sdf_file_path)
        writer.SetKekulize(False)

        for ligand in self.ligands:
            mol = ligand.mol.m  # RDKit molecule from Ligand

            properties = ligand.properties
            for prop_name, prop_value in properties.items():
                mol.SetProp(prop_name, str(prop_value))

            writer.write(mol)
        writer.close()

        return sdf_file_path

    def _to_dataframe(self):
        data = []
        for idx, ligand in enumerate(self.ligands):
            mol_props = ligand.properties

            energy_score = float(mol_props.get("Binding Energy", "0.0"))
            rscore = float(mol_props.get("Pose Score", "0.0"))

            data.append(
                {
                    "Ligand Pose Rank ID": idx + 1,
                    "Pose Score": round(rscore, 3),
                    "Binding Energy": round(energy_score, 3),
                    "Path To Docked Pose": ligand.file_path,
                }
            )
        df = pd.DataFrame(data)
        # Sort by Ranking Score descending
        df = df.sort_values(by="Pose Score", ascending=False).reset_index(drop=True)
        return df.style.format(precision=3)

    def _repr_html_(self):
        df = self._to_dataframe()
        return df._repr_html_()

    def __str__(self):
        return (
            f"DockingResult:\n  Number of Ligands: {len(self.ligands)}\n"
            f"  SMILES: {self.smiles if self.smiles else 'Not provided'}\n"
            f"  File Path: {self.file_path if self.file_path else 'Not provided'}"
        )

    def __repr__(self):
        return self.__str__()

    @jupyter_visualization
    def visualize(self, crystal_ligand_path=None, crystal_ligand_format=None):
        """
        Visualizes the docking result.
        """
        visualization_format = "sdf"
        crystal_data = None
        if crystal_ligand_format and crystal_ligand_path:
            crystal_data = {
                "raw": str(crystal_ligand_path),
                "format": crystal_ligand_format,
            }

        return DockingViewer().render_with_seperate_crystal(
            protein_data=str(self.protein.file_path),
            protein_format=self.protein.block_type,
            ligands_data=[str(self.file_path)],
            ligand_format=visualization_format,
            crystal_data=crystal_data,
        )

    def analyze(self, index: Optional[int] = None):
        """
        Analyzes the docking results and generates a summary report.
        If index is None, all ligands are analyzed.
        If index is an integer, only that specific ligand is analyzed.
        """
        if not self.ligands:
            raise ValueError("No ligands found to analyze.")

        protein = deepcopy(self.protein)
        if not protein:
            raise ValueError("No protein found to analyze.")

        fp = plf.Fingerprint()
        with tempfile.TemporaryDirectory() as temp_dir:
            protein_file = os.path.join(temp_dir, f"{protein.name}.pdb")

            protein.write_to_file(protein_file)
            protein.file_path = protein_file

            sdf_file_path = self._to_sdf(
                sdf_file_path=os.path.join(
                    temp_dir, f"{protein.name}_docking_result.sdf"
                )
            )

            v = VdWContact()
            v.vdwradii["Fe"] = 2.0
            v.vdwradii["H"] = 1.05
            v.vdwradii["O"] = 1.48

            rdkit_prot = Chem.MolFromPDBFile(protein_file, False, False)
            protein_mol = plf.Molecule(rdkit_prot)
            pose_iterable = plf.sdf_supplier(str(sdf_file_path))
            sdf_supp = Chem.SDMolSupplier(str(sdf_file_path), sanitize=False)
            pose_iterable._suppl = sdf_supp

            if index is not None:
                if index < 0 or index >= len(sdf_supp):
                    raise IndexError("Ligand index out of range.")

                single_ligand_iterable = pose_iterable[index]
                fp.run_from_iterable([single_ligand_iterable], protein_mol)

                result = fp.plot_lignetwork(single_ligand_iterable)
            else:
                fp.run_from_iterable(pose_iterable, protein_mol)
                fp.plot_barcode(xlabel="Pose")

                result = fp.to_dataframe(index_col="Pose")

        return result


class PocketFinderReport:
    def __init__(self, protein, csv_file_path=""):
        self.protein = protein
        self.file_path = csv_file_path
        self.pockets = []

    def add_pocket(self, pocket):
        self.pockets.append(pocket)

    def _to_dataframe(self):
        data = []
        for idx, pocket in enumerate(self.pockets):
            props = pocket.props
            if props:
                data.append(
                    {
                        "Pocket ID": idx + 1,
                        "Color": pocket.color,
                        "Drugability Score": props.get("drugability_score", 0),
                        "Volume": props.get("volume", 0),
                        "Total SASA": props.get("total_SASA", 0),
                        "Polar SASA": props.get("polar_SASA", 0),
                        "Polar/Apolar SASA Ratio": props.get(
                            "polar_apolar_SASA_ratio", 0
                        ),
                        "Hydrophobicity": props.get("hydrophobicity", 0),
                        "Polarity": props.get("polarity", 0),
                    }
                )

        df = pd.DataFrame(data)
        # Sort by Ranking Score descending
        df = df.sort_values(by="Drugability Score", ascending=False).reset_index(
            drop=True
        )
        return df

    def _repr_html_(self):
        df = self._to_dataframe()
        return df.style.format(precision=3)._repr_html_()

    def save_props(self):
        df = self._to_dataframe()
        df.to_csv(self.file_path, index=False)


class MolPropsReport:
    def __init__(self, results):
        self.results = results

    def _to_dataframe(self):
        smiles = []

        data = {k: [] for k in self.results[0].keys() if k != "smiles"}

        for result in self.results:
            smiles.append(result["smiles"])

            for k in data:
                data[k].append(result.get(k))

        data = {
            "SMILES": smiles,
            **data,
        }

        df = pd.DataFrame(data).style.set_properties().format(precision=3)
        return df

    def _repr_html_(self):
        df = self._to_dataframe()
        return df._repr_html_()


class ProtonationReport:
    def __init__(self, results):
        self.results, self.html_metadata = self.split_results(results)

    def split_results(self, results):
        list_without_html = []
        smiles_to_html_dict = {}

        for item in results:
            new_item = copy.deepcopy(item)
            if "html_metadata" in new_item["protonation"]:
                del new_item["protonation"]["html_metadata"]
            list_without_html.append(new_item)

            smiles = item["smiles"]
            html_meta = item["protonation"].get("html_metadata")
            smiles_to_html_dict[smiles] = html_meta

        return list_without_html, smiles_to_html_dict

    def show_plots(self):
        for smi, html_meta in self.html_metadata.items():
            centered_html = f"<center><h2>{smi}</h2>"
            display(HTML(centered_html))
            self.plot_concentration_curves(html_meta, plot=True)

    def smiles_to_img_html(self, smiles):
        mol = Chem.MolFromSmiles(smiles)
        AllChem.Compute2DCoords(mol)
        mol.GetConformer().SetId(1)

        return Draw.MolToSVG(mol, width=200, height=150).replace("\n", "")

    def plot_concentration_curves(self, html_meta, plot=False):
        fractions, smiles_list, concentration_values, pH_range, pH = html_meta
        fractions = np.transpose(np.array(fractions))

        fig = sp.make_subplots(rows=1, cols=1)
        for i, fraction in enumerate(fractions):
            fig.add_trace(
                go.Scatter(
                    x=pH_range,
                    y=fraction * 100,
                    mode="lines",
                    showlegend=True,
                    name=f"Fraction of {smiles_list[i]} at pH={pH:.1f} is {round(concentration_values[i], 2)}(%)",
                ),
                row=1,
                col=1,
            )
        fig.update_layout(
            xaxis=dict(title="pH", range=[0, 14]),
            yaxis=dict(title="Fraction (%)", range=[0, 100]),
            hovermode="closest",
        )
        if plot:
            fig.show()

    def _to_dataframe(self):
        data_dict = {}

        for result in self.results:
            smiles = result["smiles"]
            smiles_list = result["protonation"]["smiles_list"]
            concentration_list = result["protonation"]["concentration_list"]
            rounded_concentration = [round(value, 2) for value in concentration_list]
            data_dict[smiles] = pd.DataFrame(
                {
                    "protonated SMILES": smiles_list,
                    "Concentration %": rounded_concentration,
                }
            )

        df = pd.concat(data_dict.values(), keys=data_dict.keys(), names=["SMILES"])
        df["Molecule Image"] = df["protonated SMILES"].apply(self.smiles_to_img_html)
        return df

    def _repr_html_(self):
        df = self._to_dataframe()
        return df.to_html(escape=False)


class PainsReport:
    def __init__(self, results):
        self.results = results

    def get_html_of_molecule(self, result):
        molecule = Chem.MolFromSmiles(result["smiles"])
        all_matches = []
        if result["PAINS"] is not None:
            for smarts in result["PAINS"]:
                atom_matches = molecule.GetSubstructMatches(Chem.MolFromSmarts(smarts))
                all_matches.extend(atom_matches[0])

        Draw.DrawingOptions.atomHighlightsAreCircles = True
        Draw.DrawingOptions.atomHighlightColors = {
            i: (1, 0, 0) for i in set(all_matches)
        }

        img = Draw.MolToImage(molecule, size=(200, 100), highlightAtoms=all_matches)

        buffer = BytesIO()
        img.save(buffer, format="PNG")
        img_str = base64.b64encode(buffer.getvalue()).decode("utf-8")
        html = '<img src="data:image/png;base64,{0}">'.format(img_str)
        return html

    def _to_dataframe(self):
        data_dict = {}

        all_smiles_list = []
        all_molecule_html_list = []
        PAINS_pattern_list = []

        for i, result in enumerate(self.results):
            all_smiles_list.append(result["smiles"])
            all_molecule_html_list.append(self.get_html_of_molecule(result))
            PAINS_pattern_list.append(result["PAINS"])

        return pd.DataFrame.from_dict(
            {
                "SMILES": all_smiles_list,
                "Molecule Image": all_molecule_html_list,
                "SMARTS patterns of PAINS": PAINS_pattern_list,
            }
        )

    def _repr_html_(self):
        df = self._to_dataframe()
        return df.to_html(escape=False)

    def __str__(self):
        return f"DockingReport:\n  Number of DockingResults: {len(self.results)}"

    def __repr__(self):
        return self.__str__()


class DockingReport:
    def __init__(self, results: List[DockingResult], pocket_data):
        self.results = results
        self.pocket_data = pocket_data

    def _to_dataframe(self, include_props=None):
        data = []
        for result in self.results:
            property_dict = {
                "Image": None,
                "SMILES": result.smiles,
                "Ranking Score": None,
                "Binding Energy": None,
                "Path To Docked Pose": None,
            }

            if result.top_ligand and result.successful:
                ligand = result.top_ligand

                mol_props = ligand.properties

                energy_score = float(mol_props.get("Binding Energy", "0.0"))
                ranking_score = float(mol_props.get("Ranking Score", "0.0"))
                property_dict["Image"] = ligand.mol._draw()
                property_dict["SMILES"] = result.smiles
                property_dict["Ranking Score"] = round(ranking_score, 3)
                property_dict["Binding Energy"] = round(energy_score, 3)
                property_dict["Path To Docked Pose"] = ligand.file_path

                if include_props:
                    for prop in mol_props:
                        if "smiles" not in prop:
                            p = mol_props.get(prop, None)
                            property_dict[prop] = p
            data.append(property_dict)

        df = pd.DataFrame(data)
        df = df.sort_values(by="Ranking Score", ascending=False).reset_index(drop=True)
        return df

    def _repr_html_(self):
        df = self._to_dataframe().style.format(precision=3)
        return df._repr_html_()

    def generate_custom_report(self, include_props=False):
        """
        Generates a custom report as an HTML representation of a styled DataFrame.

        Args:
            include_props (list, optional): A list of properties to include in the report.
                                            If None, all properties are included.

        Returns:
            str: An HTML string representation of the styled DataFrame with specified precision.
        """
        df = self._to_dataframe(include_props).style.format(precision=3)
        return HTML(df._repr_html_())

    def __str__(self):
        return f"DockingReport:\n  Number of DockingResults: {len(self.results)}"

    def __repr__(self):
        return self.__str__()

    def save(self, save_dir=None, safe=True):
        """
        Writes the top ligands from each DockingResult to an SDF file with properties.
        """
        top_ligands = []
        for result in self.results:
            if result.top_ligand:
                top_ligands.append(result.top_ligand)

        if not top_ligands:
            return None

        if not save_dir:
            save_dir_path = (
                Path(os.getenv("END_USER_HOME", "."))
                / f"docking_report_{datetime.now().strftime('%m-%d-%Y|%H:%M:%S')}"
            )
        else:
            save_dir_path = (
                Path(save_dir)
                / f"docking_report_{datetime.now().strftime('%m-%d-%Y|%H:%M:%S')}"
            )

        save_dir_path.mkdir(parents=True, exist_ok=True)
        sdf_file_path = save_dir_path / "docking_report_top_ligands.sdf"
        if safe and sdf_file_path.exists():
            move_file_with_extension(str(sdf_file_path), "sdf")
        else:
            remove_file(str(sdf_file_path))

        writer = Chem.SDWriter(str(sdf_file_path))
        writer.SetKekulize(False)

        for ligand in top_ligands:
            mol = ligand.mol.m  # RDKit molecule

            properties = ligand.properties
            existing_properties = ligand.mol.m.GetPropsAsDict()
            if ligand.name:
                mol.SetProp("_Name", ligand.name)
            if ligand.mol.smiles:
                mol.SetProp("_SMILES", ligand.mol.smiles)

            for prop_name, prop_value in existing_properties.items():
                mol.SetProp(prop_name, str(prop_value))

            for prop_name, prop_value in properties.items():
                mol.SetProp(prop_name, str(prop_value))

            writer.write(mol)
        writer.close()

        try:
            self.results[0].protein.write_to_file(
                str(save_dir_path / f"{self.results[0].protein.name}.pdb")
            )
        except Exception as e:
            DEFAULT_LOGGER.log_error(f"Failed to write protein to file: {e}")

        save_bounding_box(
            self.pocket_data.box_center,
            self.pocket_data.box_size,
            output_file=str(save_dir_path / "bounding_box.pdb"),
        )
        return str(save_dir_path)

    @jupyter_visualization
    def visualize(
        self,
        protein_path=None,
        protein_format=None,
        sdf_file_path=None,
        crystal_ligand_path=None,
        crystal_ligand_format=None,
    ):
        """
        Visualizes the docking report by rendering the merged structures of
        protein and ligands.

        Args:
            protein_path (str, optional): Path to the protein file.
            protein_format (str, optional): Format of the protein file (e.g., pdb).
            sdf_file_path (str, optional): Path to the ligand file in SDF format.

        Raises:
            ValueError: If `protein_path` is provided without `protein_format`.

        Returns:
            Jupyter visualization object: Rendered 3D structure of the protein-ligand complex.
        """
        if sdf_file_path is None:
            file_dir = Path(self.save(save_dir="/tmp"))
            sdf_file_path = str(file_dir / "docking_report_top_ligands.sdf")

        if protein_path is not None and protein_format is None:
            raise ValueError(
                "Please provide the protein format along with the protein path."
            )

        if protein_path is None:
            if not self.results:
                raise ValueError(
                    "No results found to extract protein information from."
                )
            protein_path = str(self.results[0].protein.file_path)
            protein_format = self.results[0].protein.block_type

        viewer = DockingViewer()

        crystal_data = None
        if crystal_ligand_path and crystal_ligand_format:
            crystal_data = {
                "raw": str(crystal_ligand_path),
                "format": crystal_ligand_format,
            }

        html_content = viewer.render_with_seperate_crystal(
            protein_data=protein_path,
            protein_format=protein_format,
            ligands_data=[sdf_file_path],
            ligand_format="sdf",
            crystal_data=crystal_data,
        )

        return html_content
