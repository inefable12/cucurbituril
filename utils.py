import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.rdmolfiles import MolFromMolFile
from rdkit.Chem.rdchem import RWMol
import py3Dmol

def get_centroid(mol):
    conf = mol.GetConformer()
    coords = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    return coords.mean(axis=0)

def rotar_mol(mol, angulo_deg):
    mol = Chem.Mol(mol)
    conf = mol.GetConformer()
    angulo_rad = np.radians(angulo_deg)
    rot_matrix = np.array([
        [np.cos(angulo_rad), -np.sin(angulo_rad), 0],
        [np.sin(angulo_rad),  np.cos(angulo_rad), 0],
        [0, 0, 1]
    ])
    centroid = get_centroid(mol)
    for i in range(mol.GetNumAtoms()):
        pos = np.array(conf.GetAtomPosition(i)) - centroid
        new_pos = np.dot(rot_matrix, pos)
        conf.SetAtomPosition(i, new_pos)
    return mol

def crear_ch2_bridge():
    ch2 = Chem.MolFromSmiles("C")
    ch2 = Chem.AddHs(ch2)
    AllChem.EmbedMolecule(ch2)
    return ch2

def construir_cucurbituril_completo(mol_data, n_unidades=6):
    base = MolFromMolFile(mol_data, removeHs=False)
    base = Chem.AddHs(base)
    AllChem.EmbedMolecule(base)
    centroid = get_centroid(base)
    conf = base.GetConformer()
    for i in range(base.GetNumAtoms()):
        pos = np.array(conf.GetAtomPosition(i)) - centroid
        conf.SetAtomPosition(i, pos)

    mol_total = RWMol()
    offset = 0
    atom_index_list = []

    for i in range(n_unidades):
        angulo = (360 / n_unidades) * i
        unidad_rotada = rotar_mol(base, angulo)
        num_atoms = unidad_rotada.GetNumAtoms()
        mol_total.InsertMol(unidad_rotada)
        atom_index_list.append(offset)
        offset += num_atoms

        ch2 = crear_ch2_bridge()
        mol_total.InsertMol(ch2)
        offset += ch2.GetNumAtoms()

    for i in range(n_unidades):
        a1 = atom_index_list[i]
        a2 = atom_index_list[(i + 1) % n_unidades]
        mol_total.AddBond(a1, a2, Chem.rdchem.BondType.SINGLE)

    mol_final = mol_total.GetMol()
    AllChem.EmbedMolecule(mol_final, randomSeed=0xC0FFEE)
    AllChem.MMFFOptimizeMolecule(mol_final, maxIters=200)
    return mol_final

def showm(mol, style='stick'):
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, randomSeed=0xf00d)
    AllChem.MMFFOptimizeMolecule(mol, maxIters=200)
    mblock = Chem.MolToMolBlock(mol)

    view = py3Dmol.view(width=600, height=500)
    view.addModel(mblock, 'mol')
    view.setStyle({style: {}})
    view.zoomTo()
    return view
