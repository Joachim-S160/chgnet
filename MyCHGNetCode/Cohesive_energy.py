"""

The purpose of this python script is to calculate the cohesive energy of a structure using pymatgen

"""

from __future__ import annotations
from chgnet.model import StructOptimizer
from pymatgen.core import Structure, Molecule
from chgnet.model import CHGNet
import numpy as np


# load model
chgnet = CHGNet.load()

def Cohesive_energy(unit_cell):
    
    # load structure
    unit_cell_structure = Structure.from_file(unit_cell)
    # print(f"structure: {unit_cell_structure}")

    # Get single molecule
    single_molecule = unit_cell_structure.copy()
    sites = single_molecule.sites
    NN = single_molecule.get_sites_in_sphere(np.array([0,0,0]),3)
    print(NN)
    # Get the molecule using the site index
    molecule = Molecule.from_sites([single_molecule.sites[]])
    # print(f"single molecule: {single_molecule}")
    
    prediction = chgnet.predict_structure(unit_cell_structure)
    key, _unit = "energy", "eV/atom"
    print(f"CHGNet-predicted {key} [{_unit}]={prediction[key[0]]}\n")
    
Cohesive_energy("chgnet/MyCHGNetCode/cif_files_halides/TiBr4_mp-569814_computed.cif") 
    