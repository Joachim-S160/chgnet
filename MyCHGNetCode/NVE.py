from __future__ import annotations
from chgnet_wNVE.chgnet.model import MolecularDynamics
from chgnet_wNVE.chgnet.model import StructOptimizer
from pymatgen.core import Structure
from chgnet_wNVE.chgnet.model import CHGNet
from ase.io import read, write
from pymatgen.io.ase import AseAtomsAdaptor

def NVE_simulation(molecule_name: str="Al", input_file: str="final_al.cif",  runtime: int="1000", GPU="cuda:2"):
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
    structure = Structure.from_file(input_file)
    print(f"structure: {structure}")
    
    print(f"Number of atoms: {structure.num_sites}")


    print('structure and model loaded')
    
    md1 = MolecularDynamics(atoms=structure,
                            model=chgnet,
                            ensemble="nve",
                            timestep=1, # in fs (smaller timestep because CHGNet forces are conservative, but large MD timestep can make the numeric integration non conservative)
                            trajectory="mdNVE_out_" + molecule_name + ".traj",
                            logfile="mdNVE_out_" + molecule_name + ".log",
                            loginterval=100,
                            use_device=GPU)

    print('start md1')
    md1.run(1000*runtime)  # run an MD simulation to get a liquid structure (in ps)
    print('finished md1')
    
NVE_simulation(molecule_name="Al_combined", input_file="POSCAR", runtime=1000, GPU="cuda:3")