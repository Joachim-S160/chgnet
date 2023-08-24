from __future__ import annotations
from chgnet_wNVE.chgnet.model import MolecularDynamics
from chgnet_wNVE.chgnet.model import StructOptimizer
from pymatgen.core import Structure
from chgnet_wNVE.chgnet.model import CHGNet
from ase.md.velocitydistribution import MaxwellBoltzmannDistribution, Stationary
from ase.io import read

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
    
    # Set the momenta corresponding to the given "temperature" through velocity scaling
    MaxwellBoltzmannDistribution(atoms, temperature_K=temperature,force_temp=True)
    Stationary(atoms)  # Set zero total momentum to avoid drifting

    md1 = MolecularDynamics(atoms=atoms,
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
    
NVE_simulation(molecule_name="Al_combined_relaxed", input_file="Al_combined_relaxed.cif", runtime=0.2, GPU="cuda:0")