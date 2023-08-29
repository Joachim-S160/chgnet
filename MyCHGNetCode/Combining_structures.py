

from pymatgen.core import Structure, Lattice
from pymatgen.io.vasp.outputs import Poscar
import numpy as np

# Load the structures from files or create them as needed
structure1 = Structure.from_file("Solid_Al.cif")
structure2 = Structure.from_file("Fluid_Al.cif")

borderlength = 1 #Angstrom
# Calculate the lattice parameters for the vacuum layers (twice the original x lattice parameter)
# assert structure2.lattice.abc == structure1.lattice.abc, "The two structures must have the same lattice parameters"
lattice_params_vacuum = [2 * structure1.lattice.abc[0] + borderlength, structure1.lattice.abc[1], structure1.lattice.abc[2]] # +border in Angstrom for borders
print(f"Boxsize combined structure = {lattice_params_vacuum}")

# Calculate the angles of the lattice
# assert structure1.lattice.angles == structure2.lattice.angles, "The two structures must have the same lattice angles"
angles_vacuum = structure2.lattice.angles

# Create a new lattice with vacuum layers 
lattice_with_vacuum = Lattice.from_parameters(*lattice_params_vacuum, *angles_vacuum)

# Create vacuum layers on top and bottom of structure1 and structure2
border = np.array([[borderlength/2 , 0, 0] for _ in range(structure2.cart_coords.shape[0])])
assert border.shape == structure2.cart_coords.shape, f"The border (shape = {border.shape}) vector must have the same shape as the coordinates of structure2 (shape = {structure2.cart_coords.shape}))"
structure1_with_vacuum = Structure(lattice_with_vacuum, species=structure1.species, coords=structure1.cart_coords + border, coords_are_cartesian=True)
translation = np.array([[lattice_params_vacuum[0]*(1/2) , 0, 0] for _ in range(structure2.cart_coords.shape[0])])
assert translation.shape == structure2.cart_coords.shape, "The translation vector must have the same shape as the coordinates of structure2"
structure2_with_vacuum = Structure(lattice_with_vacuum, species=structure2.species, coords=structure2.cart_coords - border + translation, coords_are_cartesian=True)


# Append atoms from structure2_with_vacuum to structure1_with_vacuum
for site in structure2_with_vacuum:
    structure1_with_vacuum.append(site.specie, site.frac_coords, coords_are_cartesian=False)

# Save the combined structure to a CIF file
structure1_with_vacuum.to(filename="solid_fluid_interface_Al.cif")
print("combined structure has been made")

# Melted to fast with the heating method: WCl6, TiI4, TiBr4,