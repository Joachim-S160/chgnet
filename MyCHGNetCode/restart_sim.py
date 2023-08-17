"""

Description:
restarting an NPT simulation from a previous NPT simulation

"""

from __future__ import annotations
from chgnet.model import MolecularDynamics
from chgnet.model import StructOptimizer
from pymatgen.core import Structure
from chgnet.model import CHGNet
from ase.io.trajectory import Trajectory


def simulation(file: str = "structure.traj", molecule_name: str = 'WCl6', GPU:str="cuda:1") -> None:
    """
    Args: 
        file: the trajectory file .traj
    Returns:
        None, stores .traj and .log file of new NPT run
    """

    # initiate the model
    chgnet = CHGNet.load()
    # get the structure from the last frame in the trajectory file
    structure = Trajectory(file)[-1]
    print('Presimulation code works, start md1')
    # start md simulation
    md = MolecularDynamics(
        atoms=structure,
        # atoms = atoms_last_frame,
        model=chgnet,
        ensemble="npt",
        # compressibility_au=2.1,
        temperature=2000,  # in K
        timestep=2,  # in fs
        trajectory="mdNPT2_out_restart_" + molecule_name + ".traj",
        logfile="mdNPT2_out_restart_" + molecule_name + ".log",
        loginterval=100,
        use_device=GPU,  # use 'cuda:2' for faster MD
        taut=0.5*100*1000  # in fs
    )
    md.run(1*500*1000)  # run a 1 ns MD simulation
    
simulation(file="mdNPT2_out_HfF4.traj", molecule_name="HfF4", GPU="cuda:1")
