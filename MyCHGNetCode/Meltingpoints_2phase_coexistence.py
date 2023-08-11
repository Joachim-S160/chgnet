"""

Description:
Aluminum melting point prediction using CHGNet
using the 2 phase coexistence method

"""

from __future__ import annotations
from chgnet.model import MolecularDynamics
from chgnet.model import StructOptimizer
from pymatgen.core import Structure
from chgnet.model import CHGNet

# get biggest box with less than 500 atoms


def Biggest_box(structure):
    """

    input: 
        pymatgen structure
    output: 
        list of boxdimensions [a, b, c] with total number of atoms > 500

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


def simulation(molecule_name, cif_file, GPU="cuda:2"):
    """

    input:
        molecule_name: name of the molecule
        cif_file: cif file of the molecule
        GPU: GPU to use
    output:
        stores trajectory and log file in .traj and .log files

    """

    # load structure
    structure = Structure.from_file(
        f"chgnet/MyCHGNetCode/cif_files_frames/HfF4/HfF4_frame_6.cif")
    print(f"structure: {structure}")

    # load model
    chgnet = CHGNet.load()

    # relax structure at 0K
    relaxer = StructOptimizer()
    relaxed_structure = relaxer.relax(structure, verbose=True)
    print(
        f"\nCHGNet took {len(relaxed_structure['trajectory'])} steps. Relaxed structure:")
    print(relaxed_structure["final_structure"])

    # create solid supercell
    Solid:structure = structure.make_supercell(Biggest_box(structure))

    # molecular dynamics simulations
    # melt at 2000K via nvt
    md1 = MolecularDynamics(atoms=Solid,
                            model=chgnet,
                            ensemble="nvt",
                            temperature=2000, # in K
                            timestep=2, # in fs
                            trajectory="mdNVT_out_" + molecule_name + ".traj",
                            logfile="mdNVT_out_" + molecule_name + ".log",
                            loginterval=100,
                            use_device=GPU)
    
    print('start md1')
    md1.run(500*1000)  # run a 1 ns MD simulation to get a liquid structure
    print('finished md1')
    Liquid:structure = md1.atoms
    
    # create solid liquid interface
    for site in Solid:
        Liquid.append(site.specie, site.frac_coords, coords_are_cartesian=False)
    combined_structure:structure = Liquid
    atoms = combined_structure.atoms()
    atoms.write("combined_structure.cif", format='cif')
    print("combined structure has been made")




