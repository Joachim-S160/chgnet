"""
Purpose: Test if a melted structure has truely melted by relaxing it and checking if the atoms don't revert to the solid structure.
Only do on a supercomputer because it takes a long time to run when using +1000 atoms.
"""

from plot import get_molecule_name
from ase.io.trajectory import Trajectory
from chgnet.model import MolecularDynamics
from chgnet.model import StructOptimizer
from pymatgen.core import Structure
from chgnet.model import CHGNet
import numpy as np

def Melted_or_not(ciffile:str = 'Al.cif'):
    
    chgnet = CHGNet.load()
    structure = Structure.from_file(ciffile)
    # Relax the structure so that the atoms are moved to positions with lower potential energy and the cell size is adjusted to the optimal size with no stresses.
    relaxer = StructOptimizer()
    relaxed_structure_dict:dict = relaxer.relax(structure, verbose=True)
    relaxed_structure_dict["final_structure"].to(filename="Test_melted" + ".cif")
    
Melted_or_not("chgnet/MyCHGNetCode/2PC_strcutures/LiCL/Liquid.cif")
    