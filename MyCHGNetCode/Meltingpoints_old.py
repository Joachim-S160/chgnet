"""

Description:
Aluminum melting point prediction using CHGNet
using the heating curve method
1 NVT run at 400K
1 NPT heating run to 2000K

"""

from __future__ import annotations
from chgnet.model import MolecularDynamics
from chgnet.model import StructOptimizer
from pymatgen.core import Structure
from chgnet.model import CHGNet


# load structure
structure = Structure.from_file(f"Al.cif")
print(f"structure: {structure}")

# load model
chgnet = CHGNet.load()

# Molecular dynamics simulations

# Make a 7x7x7 supercell structure which is 343 copies of the original structure
structure.make_supercell([7,7,7])

# # eq at 400K via nvt
# md1 = MolecularDynamics(
#     atoms=structure,
#     model=chgnet,
#     ensemble="nvt",
#     temperature=400,  # in K
#     timestep=2,  # in fstaut=500*1000000
#     trajectory="mdNVT_out.traj",
#     logfile="mdNVT_out.log",
#     loginterval=100,
#     use_device="cuda",  # use 'cuda' for faster MD
# )

# md1.run(500*20)

from ase.io.trajectory import Trajectory
# get last frame of md1
atoms_last_frame = Trajectory("mdNVT_out.traj")[-1]


md2 = MolecularDynamics(
    # atoms=md1.atoms,
    atoms = atoms_last_frame,
    model=chgnet,
    ensemble="npt",
    compressibility_au=2.1,
    temperature=2000,  # in K
    timestep=2,  # in fs
    trajectory="mdNPT2_out.traj",
    logfile="mdNPT2_out.log",
    loginterval=100,
    use_device="cuda:2",  # use 'cuda:2' for faster MD
    taut=0.5*100*1000 #in fs
)
md2.run(5*500*1000)  # run a 0.1 ps MD simulation

