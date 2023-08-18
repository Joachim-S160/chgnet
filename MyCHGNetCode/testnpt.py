"""

small npt test to check out the ideal gas contribution to the pressure from the ase code

"""

from __future__ import annotations
from chgnet.model import MolecularDynamics
from chgnet.model import StructOptimizer
from pymatgen.core import Structure
from chgnet.model import CHGNet

def simulation_0():
    

    cif_file = "WCl6.cif"

    # load structure
    structure = Structure.from_file(cif_file)
    print(f"structure: {structure}")

    # load model
    chgnet = CHGNet.load()

    # Relax the structure so that the atoms are moved to positions with lower potential energy and the cell size is adjusted to the optimal size with no stresses.
    relaxer = StructOptimizer()
    relaxed_structure_dict:dict = relaxer.relax(structure, verbose=True)
    print(
        f"\nCHGNet took {len(relaxed_structure_dict['trajectory'])} steps. Relaxed structure:")
    print(relaxed_structure_dict["final_structure"])
    relaxed_structure:Structure = relaxed_structure_dict["final_structure"]

    # Check if relaxed structure is a structure type
    assert isinstance(relaxed_structure,Structure), "Relaxed structure is not a structure type"

    # Molecular dynamics simulations
    name = "default_with_ideal_gas_contributions"
    print('Presimulation code works, start md1')
    # eq at 400K via nvt
    md1 = MolecularDynamics(
        atoms=relaxed_structure,
        model=chgnet,
        ensemble="npt",
        temperature=100,  # in K
        timestep=2,  # in fs, taut=500*1000000
        trajectory="mdNPT_out_" + name + ".traj",
        logfile="mdNPT_out_" + name + ".log",
        loginterval=100,
        use_device="cuda",  # use 'cuda' for faster MD
    )
    md1.run(500*20)  # run 20 ps


from plot import time_pressure_plot, time_temperature_plot
time_pressure_plot(["chgnet/MyCHGNetCode/data_out_test/mdNPT_out_default_with_ideal_gas_contributions.traj","chgnet/MyCHGNetCode/data_out_test/mdNPT_out_without_ideal_gas_contributions.traj"])
# time_temperature_plot(["chgnet/MyCHGNetCode/data_out_test/mdNPT_out_default_with_ideal_gas_contributions.traj","chgnet/MyCHGNetCode/data_out_test/mdNPT_out_without_ideal_gas_contributions.traj"])