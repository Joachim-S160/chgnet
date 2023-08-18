"""

Description:
Halide melting point prediction using CHGNet
using the heating curve method
1 NVT thermalisation at Tstart
1 NPT Equilibration at Tstart
1 NPT heating run to Tend
no more than 500 atoms in the box
WATCH OUT WITHOUT RELAXING THE STRUCTURE FIRST THE MD CAN EXPLODE

"""

from __future__ import annotations
from chgnet.model import MolecularDynamics
from chgnet.model import StructOptimizer
from pymatgen.core import Structure
from chgnet.model import CHGNet

# get biggest box with less than 500 atoms

def Biggest_box(structure):
    # Get number of atoms in unitcell
    noau = structure.num_sites
    print(noau)
    # Get the boxdimension
    Boxdim = int((500/noau)**(1/3))
    a, b, c = Boxdim, Boxdim, Boxdim

    noa = (a+1)*b*c*noau
    # Update boxdimensions until noa is bigger than 500 or the boxdimensions are the same as before
    while noa < 500:
        x, y, z = a, b, c
        if (a+1)*b*c*noau < 500:
            a += 1

        if a*(b+1)*c*noau < 500:
            b += 1

        if a*b*(c+1)*noau < 500:
            c += 1
        if [x, y, z] == [a, b, c]:
            break
    return [a, b, c]


def Melting_point_simulation(molecule_name, cif_file, Tstart=300, Tend=2000, GPU="cuda:2"):

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

    # Make a axbxc supercell structure which is a*b*c copies of the original structure (noau atoms) = a*b*c*noau atoms
    structure.make_supercell(Biggest_box(structure))

    print('Presimulation code works, start md1')
    # eq at 400K via nvt
    md1 = MolecularDynamics(
        atoms=relaxed_structure,
        model=chgnet,
        ensemble="nvt",
        temperature=Tstart,  # in K
        timestep=2,  # in fs, taut=500*1000000
        trajectory="mdNVT_out_" + molecule_name + ".traj",
        logfile="mdNVT_out_" + molecule_name + ".log",
        loginterval=100,
        use_device=GPU,  # use 'cuda' for faster MD
    )

    md1.run(500*10)  # run a 10 ps MD simulation
    print('md1 works, start md2')

    md2 = MolecularDynamics(
        atoms=md1.atoms,
        model=chgnet,
        ensemble="npt",
        temperature=Tstart,  # in K
        timestep=2,  # in fs, taut=500*1000000
        trajectory="mdNPT1_out_" + molecule_name + ".traj",
        logfile="mdNPT1_out_" + molecule_name + ".log",
        loginterval=100,
        use_device=GPU,  # use 'cuda' for faster MD
        taut=0.5*100*1000  # in fs
    )

    md2.run(500*10)  # run a 10 ps MD simulation
    print('md2 works, start md3')

    # from ase.io.trajectory import Trajectory
    # # get last frame of md1
    # atoms_last_frame = Trajectory("mdNVT_out.traj")[-1]

    md3 = MolecularDynamics(
        atoms=md2.atoms,
        # atoms = atoms_last_frame,
        model=chgnet,
        ensemble="npt",
        # compressibility_au=2.1,
        temperature=Tend,  # in K
        timestep=2,  # in fs
        trajectory="mdNPT2_out_" + molecule_name + ".traj",
        logfile="mdNPT2_out_" + molecule_name + ".log",
        loginterval=100,
        use_device=GPU,  # use 'cuda:2' for faster MD
        taut=0.5*100*1000  # in fs
    )
    md3.run(1*500*1000)  # run a 1 ns MD simulation
    print('md3 works, simulation finished')

Melting_point_simulation('TbCl3', 'TbCl3.cif', 300, 1500, "cuda:3")

# Ask Bowen what default pressure is used in the NPT ensemble, show him pressure plot