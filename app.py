import requests
from rdkit import Chem
from rdkit.Chem import AllChem

def cargar_glicolurilo_desde_github(url_xyz):
    response = requests.get(url_xyz)
    if response.status_code != 200:
        raise Exception("No se pudo descargar el archivo .xyz desde GitHub.")
    
    xyz_block = response.text
    mol = Chem.MolFromXYZBlock(xyz_block)
    if mol is None:
        raise Exception("No se pudo convertir el archivo .xyz a una molécula válida con RDKit.")
    
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    return mol
