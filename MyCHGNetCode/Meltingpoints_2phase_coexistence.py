"""

Description:
Aluminum melting point prediction using CHGNet
and the 2 phase coexistence method

"""

from __future__ import annotations
from chgnet.model import MolecularDynamics
from chgnet.model import StructOptimizer
from pymatgen.core import Structure
from chgnet.model import CHGNet
from ase.io import read, write
from pymatgen.io.ase import AseAtomsAdaptor

def Biggest_box(structure):
    """

    input: 
        pymatgen structure
    output: 
        list of boxdimensions [a, b, c] with total number of atoms < 500

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

    # relax structure at 0K
    relaxer = StructOptimizer()
    relaxed_structure_dict:dict = relaxer.relax(structure, verbose=True)
    print(
        f"\nCHGNet took {len(relaxed_structure_dict['trajectory'])} steps. Relaxed structure:")
    print(relaxed_structure_dict["final_structure"])
    relaxed_structure:Structure = relaxed_structure_dict["final_structure"]
    
    # Check if relaxed structure is a structure type
    assert isinstance(relaxed_structure,Structure), "Relaxed structure is not a structure type"
    
    # create solid supercell
    relaxed_structure.make_supercell(Biggest_box(relaxed_structure))
    Solid = relaxed_structure

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
    Solid.to(filename="Solid_" + molecule_name + ".cif") #Solid:Structure, so .to() works and .write() doesn't
    Liquid.write(filename="Liquid_" + molecule_name + ".cif", format='cif') #Liquid:md.atoms,
    print("Solid and liquid cif files have been made")

def NVE_simulation(input_file: str="final_al.cif", molecule_name: str="Al", runtime: int="1000", GPU="cuda:2"):
    """
    Args:
        input_file: input file of the molecule
        molecule_name: name of the molecule
        runtime: runtime of the NVE simulation in ps
        GPU: GPU to use
    Returns:
        None. Stores trajectory and log files in .traj and .log files
    """
    
    # load structure
    Atoms = read(input_file)
    structure = AseAtomsAdaptor.get_structure(Atoms)
    # print(f"structure: {structure}")

    # load model
    chgnet = CHGNet.load()
    print('structure and model loaded')
    
    md1 = MolecularDynamics(atoms=structure,
                            model=chgnet,
                            ensemble="nve",
                            timestep=2, # in fs
                            trajectory="mdNVE_out_" + molecule_name + ".traj",
                            logfile="mdNVE_out_" + molecule_name + ".log",
                            loginterval=100,
                            use_device=GPU)

    print('start md1')
    md1.run(500*runtime)  # run an MD simulation to get a liquid structure (in ps)
    print('finished md1')
    


# simulation(molecule_name="Al", cif_file="Al.cif", temperature_fluid=1800, NVT_runtime=1000, GPU="cuda:2")
NVE_simulation(molecule_name="Al_combined", input_file="Al_final.cfg", runtime=0.1, GPU="cuda:2")

