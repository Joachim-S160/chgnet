"""

Description:
Melting point prediction using CHGNet
and the 2 phase coexistence method
WATCH OUT WITHOUT RELAXING THE STRUCTURE FIRST THE MD CAN EXPLODE

"""

from __future__ import annotations
from chgnet.model import MolecularDynamics
from chgnet.model import StructOptimizer
from pymatgen.core import Structure
from chgnet.model import CHGNet
from ase.io import read, write
from pymatgen.io.ase import AseAtomsAdaptor

def Biggest_elongated_box(structure=None, a=5, b=1, c=1, x_factor=5, Max_atoms=600):
    """
    Args:   
        structure: pymatgen structure 
        a,b,c: initial box dimensions or lattice parameters,
        x_factor: ratio between a and max(b,c) with a,b,c the lattice parameters
    Returns: 
        dimensions of the biggest possible supercell box  
    
    Note this function prefers the stretching of the box in the x direction
    Has extra functionality built in so that b and c lattice parameters are always bigger than 10 Angstrom
    """
    # Get number of atoms in unitcell

    noau = structure.num_sites
    A_lattice, B_lattice, C_lattice = structure.lattice.abc

    print(f"Number of atoms in unitcell: ", noau)

    noa = a * b * c * noau
    assert noa != 0, "No atoms in unitcell"
    assert noa < Max_atoms, f"Box is already bigger than {Max_atoms} atoms"
    # Update boxdimensions until noa is bigger than 500 or the boxdimensions are the same as before
    while noa < Max_atoms:
        x, y, z = a, b, c
        if (a + 1) * b * c * noau < Max_atoms:
            a += 1

        if a * (b + 1) * c * noau < Max_atoms and a >= x_factor * max(b + 1, c):
            b += 1

        if a * b * (c + 1) * noau < Max_atoms and a >= x_factor * max(b, c + 1):
            c += 1
        if [x, y, z] == [a, b, c]:
            break
        noa = a * b * c * noau
    print("Box dimensions: ", a, b, c)
    print("Total number of atoms in the box: ", a*b*c*noau)
    print(f"Boxsize supercell structure = {a * A_lattice, b * B_lattice, c * C_lattice} [Angstrom]")
    assert b * B_lattice >= 10, f"b lattice parameter is smaller than 10 Angstrom, b lattice parameter is {b * B_lattice}"
    assert c * C_lattice >= 10, f"c lattice parameter is smaller than 10 Angstrom, c lattice parameter is {c * C_lattice}"
    assert a *A_lattice * 5 >= b * B_lattice, f"Box is not elongated enough, a lattice parameter is {a * A_lattice} and b lattice parameter is {b * B_lattice}"
    assert a *A_lattice * 5 >= c * C_lattice, f"Box is not elongated enough, a lattice parameter is {a * A_lattice} and c lattice parameter is {c * C_lattice}"

    return [a, b, c]

def simulation(molecule_name:str="Al", cif_file:str="Al.cif", temperature_fluid:int=2000, NVT_runtime:int="1000", GPU="cuda:2"):
    """

    input:
        molecule_name: name of the molecule
        cif_file: cif file of the molecule
        temperature_fluid: temperature of the fluid in K
        NVT_runtime: runtime of the NVT simulation in ps
        GPU: GPU to use
    output:
        stores trajectory and log files in .traj and .log files
        stores cif files of the solid @ 0K and fluid @ temperature_fluid K, but NOT the combined structure yet

    """

    # load structure
    structure = Structure.from_file(cif_file)
    print(f"structure: {structure}")

    # load model
    chgnet = CHGNet.load()

    # Relax the structure so that the atoms are moved to positions with lower potential energy and the cell size is adjusted to the optimal size with no stresses.
    relaxer = StructOptimizer(use_device=GPU)
    relaxed_structure_dict:dict = relaxer.relax(structure, verbose=True)
    print(
        f"\nCHGNet took {len(relaxed_structure_dict['trajectory'])} steps. Relaxed structure:")
    print(relaxed_structure_dict["final_structure"])
    relaxed_structure:Structure = relaxed_structure_dict["final_structure"]
    
    # Check if relaxed structure is a structure type
    assert isinstance(relaxed_structure,Structure), "Relaxed structure is not a structure type"
    
    # create solid supercell
    relaxed_structure.make_supercell(Biggest_elongated_box(relaxed_structure))
    Solid = relaxed_structure

    # Test if the solid has at least 100 atoms
    assert Solid.num_sites >= 100, "Solid has less than 100 atoms, too much finite size error will occur"

    # create cif file of solid
    Solid.to(filename="Solid_" + molecule_name + ".cif") #Solid:Structure, so .to() works and .write() doesn't

    # molecular dynamics simulations
    # melt at 2000K via nvt
    md1 = MolecularDynamics(atoms=Solid,
                            model=chgnet,
                            ensemble="nvt",
                            temperature=temperature_fluid, # in K
                            timestep=2, # in fs
                            trajectory="mdNVT_out_" + molecule_name + ".traj",
                            logfile="mdNVT_out_" + molecule_name + ".log",
                            loginterval=100,
                            use_device=GPU)

    print('start md1')
    md1.run(500*NVT_runtime)  # run an MD simulation to get a liquid structure (in ps)
    print('finished md1')
    Liquid = md1.atoms #ase.Atoms object

    # create cif files of solid and liquid
    
    Liquid.write(filename="Liquid_" + molecule_name + ".cif", format='cif') #Liquid:md.atoms,
    print("Solid and liquid cif files have been made")

simulation(molecule_name="LiF", cif_file="LiF.cif", temperature_fluid=1800, NVT_runtime=30, GPU="cuda:2")


