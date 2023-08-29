from __future__ import annotations
from chgnet.model import MolecularDynamics
from chgnet.model import StructOptimizer
from pymatgen.core import Structure, Lattice
from pymatgen.io.ase import AseAtomsAdaptor
from chgnet.model import CHGNet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.io import read
import numpy as np
import os

def create_directory(directory_path):
    """
    Args:
        directory_path: path to the directory to be created
    Returns:
        None. Creates directory if it does not exist yet    
    """
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
        print(f"Directory '{directory_path}' created.")
    else:
        print(f"Directory '{directory_path}' already exists.")


# expand_volume(Structure.from_file("chgnet/MyCHGNetCode/2PC_strcutures/LiCL/LiCl_combined_relaxed.cif"), "LiCl", factor=27)

def NVE_simulation(molecule_name: str="Al", input_file: str="final_al.cif",  runtime: int=1000, temperature:int=1000, GPU="cuda:2"):
    """
    Args:
        input_file: input file of the molecule
        molecule_name: name of the molecule
        runtime: runtime of the NVE simulation in ps
        GPU: GPU to use
    Returns:
        None. Stores trajectory and log files in .traj and .log files
    """
  
    # load model
    chgnet = CHGNet.load()    
    # load structure
    atoms = read(input_file) #ASE structure
    print('structure and model loaded')

    #######################Uncomment to have expansion and relaxation in one script#######################
    
    # Expand the structure:
    structure = AseAtomsAdaptor.get_structure(atoms)
    print('Expand structure by a factor of 1.03')
    a = 0.03**(1/3)
    structure.apply_strain(np.array([a,a,a]))   
    print('Structure expansion finished') 
    
    # Check with ovito if the structure is expanded correctly
    create_directory(f"{os.getcwd()}/expanded_structure")
    structure.to(filename=f"{os.getcwd()}/expanded_structure/{molecule_name}_{1.03}_percent.cif")

    # Perform ionic relaxation (no Volume change)
    # = Relax the structure so that the atoms are moved to positions with lower potential energy but the cell size is NOT adjusted to the optimal size with no stresses.
    relaxer = StructOptimizer(use_device=GPU)
    relaxed_structure_dict: dict = relaxer.relax(
        structure, fmax=0.1, relax_cell=False)
    relaxed_structure: Structure = relaxed_structure_dict["final_structure"]
    print("Relaxation finished")
    
    # Check with ovito if the structure is relaxed correctly
    relaxed_structure.to(filename=molecule_name + "_relaxed" + ".cif")
    print('relaxed structure saved')

    #####################################################################################################
    
    atoms1 = AseAtomsAdaptor.get_atoms(relaxed_structure)
    # Set the momenta corresponding to the given "temperature" through velocity scaling
    MaxwellBoltzmannDistribution(atoms1, temperature_K=temperature,force_temp=True)
    Stationary(atoms1)  # Set zero total momentum to avoid drifting

    md1 = MolecularDynamics(atoms=atoms1,
                            model=chgnet,
                            ensemble="nve",
                            timestep=1, # in fs (smaller timestep because CHGNet forces are conservative, but large MD timestep can make the numeric integration non conservative)
                            trajectory="mdNVE_out_" + molecule_name + ".traj",
                            logfile="mdNVE_out_" + molecule_name + ".log",
                            loginterval=100,
                            use_device=GPU)

    print('start md1')
    md1.run(1000*runtime)  # run an MD simulation (in ps)
    print('finished md1')
    
NVE_simulation(molecule_name="LiF_combined_relaxed_expanded", input_file="LiF_junction_relaxed.cif", runtime=200, GPU="cuda:3")