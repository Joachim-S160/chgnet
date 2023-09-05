"""

The purpose of this script is to combine two structures into one structure using pymatgen

"""

from pymatgen.core import Structure, Lattice
from pymatgen.io.vasp.outputs import Poscar
import numpy as np

def combined_structure(Solid_cif:str="Solid_Al.cif", Fluid_cif:str="Fluid_Al.cif", name:str="solid_fluid_interface_LiS.cif", separation_vacuum_length:float=2):
    """
    Args:
        Solid_cif: path to the cif file of the solid structure
        Fluid_cif: path to the cif file of the fluid structure
        name: name of the combined structure
        separation_vacuum_length: length of the vacuum layer between the two structures at the interfaces [Angstrom]
    Returns:    
        None. Stores the combined structure in a cif file
    Note: 
        The two structures must have the same lattice parameters and angles, which after an NVT simulation is the case
    """

    # Load the structures from files or create them as needed
    Solid = Structure.from_file(Solid_cif)
    Fluid = Structure.from_file(Fluid_cif)

    # Calculate the angles of the lattice
    # assert Solid.lattice.angles == Fluid.lattice.angles, "The two structures must have the same lattice angles"
    angles_vacuum = Fluid.lattice.angles

    # Calculate the lattice parameters for the vacuum layers (twice the original x lattice parameter)
    # assert Fluid.lattice.abc == Solid.lattice.abc, "The two structures must have the same lattice parameters"
    lattice_params_vacuum = [2 * Fluid.lattice.abc[0] + 2 * separation_vacuum_length, 
                             Fluid.lattice.abc[1], Fluid.lattice.abc[2]] 
    print(f"Boxsize combined structure = {lattice_params_vacuum}")
    print(f"Boxsize solid structure = {Solid.lattice.abc}")
    print(f"Boxsize fluid structure = {Fluid.lattice.abc}")
    
    # Create a new lattice with nothing in it
    lattice_with_vacuum = Lattice.from_parameters(*lattice_params_vacuum, *angles_vacuum)
    
    # Add the atoms from the solid structure to a copy of the new lattice
    Vacuum_with_solid = Structure(lattice_with_vacuum, species=Solid.species, coords=Solid.cart_coords, coords_are_cartesian=True)
    
    # Create translation vectors for the fluid structure
    translation = np.array([[0.5 , 0, 0] for _ in range(Fluid.cart_coords.shape[0])])
    # Create new coordinates for the fluid structure
    new_coords = np.zeros_like(Fluid.frac_coords)
    new_coords[:,0] += translation[:,0]
    # Adding vacuum by scaling
    new_coords[:,0] += Fluid.frac_coords[:,0]/2 * (2 * Fluid.lattice.abc[0] / lattice_params_vacuum[0])
    new_coords[:,1] += Fluid.frac_coords[:,1]
    new_coords[:,2] += Fluid.frac_coords[:,2]
    # Add the atoms from the fluid structure to a copy of the new lattice
    Vacuum_with_fluid = Structure(lattice_with_vacuum, species=Fluid.species, 
                                  coords=new_coords, coords_are_cartesian=False)
    
    # Append atoms from Vacuum_with_fluid to Vacuum_with_solid
    for site in Vacuum_with_fluid:
        Vacuum_with_solid.append(site.specie, site.coords, coords_are_cartesian=True)
        
    # Save the combined structure to a CIF file
    Vacuum_with_solid.to(filename=name)
    print("combined structure has been made")

combined_structure(Solid_cif="Solid_LiCl.cif", Fluid_cif="Fluid_LiCl.cif", name = "solid_fluid_interface_LiCl.cif", separation_vacuum_length=2)