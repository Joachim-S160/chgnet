"""

Description:
Halide melting point prediction using CHGNet
using the dicrete temperature density method
1 NPT Equilibration at different temerpatures
no more than 500 atoms in the box

"""

from __future__ import annotations
from chgnet.model import MolecularDynamics
from chgnet.model import StructOptimizer
from pymatgen.core import Structure
from chgnet.model import CHGNet

# get biggest box with less than 500 atoms


def Biggest_box(structure):
    """
    Args: structure
    Returns: dimensions of the biggest possible supercell box    
    """
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


def Melting_point_simulation(molecule_name, cif_file, Tstart=100, Tend=1000, GPU="cuda:2"):
    """
    Args:   molecule_name: name of the molecule
            cif_file: cif file of the molecule
            Tstart: start temperature
            Tend: end temperature
            GPU: GPU to use
    Returns: None, stores trajectory and log file
    """

    # load structure
    structure = Structure.from_file(cif_file)
    print(f"structure: {structure}")

    # load model
    chgnet = CHGNet.load()

    # Molecular dynamics simulations

    # Make a axbxc supercell structure which is a*b*c copies of the original structure (noau atoms) = a*b*c*noau atoms
    structure.make_supercell(Biggest_box(structure))

    print('Presimulation code works, start md1')
    for temp in range(Tstart, Tend, 100):
        print(temp)
        md1 = MolecularDynamics(
            atoms=structure,
            model=chgnet,
            ensemble="npt",
            temperature=temp,  # in K
            timestep=2,  # in fs, taut=500*1000000
            trajectory="mdNPT_out_" + molecule_name + \
            "_" + str(temp) + ".traj",
            logfile="mdNPT_out_" + molecule_name + "_" + str(temp) + ".log",
            loginterval=100,
            use_device=GPU,  # use 'cuda' for faster MD
        )
        md1.run(500*40)  # run 40 ps
        print('md ' + str(temp) + ' done')

Melting_point_simulation("WCl6", "WCl6.cif", Tstart=100, Tend=1000, GPU="cuda:2")