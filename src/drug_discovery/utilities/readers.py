import os
import numpy as np

from rdkit import Chem
from .utils import Logger
from ..structures.internal_structures import FileFormat


N_PT_ELEMENTS = 118
PERIODIC_TABLE = {
    'id': [],
    'symbol': [],
    'name': [],
    'weight': [],
    'rvdw': [],
    'rcov': []
}

def init_periodic_table():
    logger = Logger("INFO", os.getenv("LOG_BIOSIM_CLIENT"))
    pt = Chem.GetPeriodicTable()

    for id in range(1, N_PT_ELEMENTS + 1):
        try:
            symbol = pt.GetElementSymbol(id)
            name = pt.GetElementName(id)
            weight = pt.GetAtomicWeight(id)
            rvdw = pt.GetRvdw(id)
            rcov = pt.GetRcovalent(id)

            assert isinstance(symbol, str), "Symbol must be a string"
            assert isinstance(name, str), "Name must be a string"
            assert isinstance(weight, float), "Weight must be a float"
            assert isinstance(rvdw, float), "rvdw must be a float"
            assert isinstance(rcov, float), "rcov must be a float"

            PERIODIC_TABLE['id'].append(id)
            PERIODIC_TABLE['symbol'].append(symbol)
            PERIODIC_TABLE['name'].append(name)
            PERIODIC_TABLE['weight'].append(weight)
            PERIODIC_TABLE['rvdw'].append(rvdw)
            PERIODIC_TABLE['rcov'].append(rcov)
        except Exception as e:
            logger.log_error(f"Failed to load Periodic Table on atom with id: {id}. Error: {str(e)}")


def from_pdbqt_atom_type(atom_type: str):
    if atom_type in ["A", "Z", "G", "GA", "J", "Q"]:
        return "C"
    elif atom_type in ["HD", "HS"]:
        return "H"
    elif atom_type in ["NA", "NS"]:
        return "N"
    elif atom_type in ["OA", "OS"]:
        return "O"
    elif atom_type == "SA":
        return "S"
    else:
        return ""
    

def read_pdb_pdbqt_block(block_type, block):
    lines = block.split("\n")
    uppercase_atom_symbols = [symbol.upper() for symbol in PERIODIC_TABLE['symbol']]

    if block_type == FileFormat.PDB:
        compnd_line_id = next((i for i, line in enumerate(lines) if line.startswith("COMPND")), None)
        name = lines[compnd_line_id][7:].strip() if compnd_line_id is not None else None
    else:
        name_line_id = next((i for i, line in enumerate(lines) if line.startswith("REMARK  Name")), None)
        name = lines[name_line_id][15:].strip() if name_line_id is not None else None

    atom_line_ids = [i for i, line in enumerate(lines) if line.startswith("ATOM")]

    atom_count = len(atom_line_ids)

    atom_types = [None] * atom_count
    coordinates = np.empty((3, atom_count), dtype=np.float32)

    atom_lines = [lines[i] for i in atom_line_ids]
    for i, line in enumerate(atom_lines):
        coordinates[0, i] = float(line[30:38])
        coordinates[1, i] = float(line[38:46])
        coordinates[2, i] = float(line[46:54])

        atom_type = line[70:73].strip()
        if atom_type not in uppercase_atom_symbols:
            atom_type = from_pdbqt_atom_type(atom_type)

        if atom_type == '':
            atom_type = line[76:79].strip()
            if atom_type not in uppercase_atom_symbols:
                atom_type = from_pdbqt_atom_type(atom_type)

        atom_types[i] = atom_type

    return name, atom_types, coordinates


def read_mol2_block(block: str):
    mol = Chem.rdmolfiles.MolFromMol2Block(block, removeHs = False, sanitize = False)
    name = mol.GetProp('_Name')
    
    atoms = mol.GetAtoms()
    atom_types = np.empty(len(atoms), dtype=str)
    coordinates = np.empty((3, len(atoms)), dtype=np.float32)

    for i, atom in enumerate(atoms):
        atom_types[i] =  atom.GetSymbol()
        positions = mol.GetConformer(0).GetAtomPosition(i)
        coordinates[:, i] = np.array([positions.x, positions.y, positions.z], dtype=np.float32)

    return name, atom_types, np.transpose(coordinates)


def read_xyz_block(block: str):
    lines = block.split("\n")

    atom_count = int(lines[0])
    name = lines[1].strip()

    atom_types = np.empty(atom_count, dtype=str)
    coordinates = np.empty((3, atom_count), dtype=np.float32)

    atom_lines = lines[2:2+atom_count]
    for i, line in enumerate(atom_lines):
        tokens = line.split()
        atom_types[i] = tokens[0]
        coordinates[0, i] = float(tokens[1])
        coordinates[1, i] = float(tokens[2])
        coordinates[2, i] = float(tokens[3])

    return name, atom_types, coordinates

def read_block(block_type, block_content):
    if block_type == FileFormat.MOL2:
        return read_mol2_block(block_content)
    elif block_type in [FileFormat.PDB, FileFormat.PDBQT]:
        return read_pdb_pdbqt_block(block_type, block_content)
    elif block_type == FileFormat.XYZ:
        return read_xyz_block(block_content)
    else:
        raise Exception(f"Invalid file format {block_type}")