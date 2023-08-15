from pymatgen.core import Structure, Lattice
from pymatgen.io.vasp.outputs import Poscar
import numpy as np

# Load your structures from files or create them as needed
structure1 = Structure.from_file("TiBr4.cif")
structure2 = Structure.from_file("Tl2Cl3FAKE.cif")

# Calculate the lattice parameters for the combined structure (same x lattice parameter, 4 times larger in y and z)
lattice_params_vacuum = [structure1.lattice.abc[0], structure1.lattice.abc[1], structure1.lattice.abc[2]]
print(f"Boxsize combined structure = {lattice_params_vacuum}")

# Calculate the angles of the lattice
angles_vacuum = structure1.lattice.angles

# Create a new lattice with vacuum layers on top and bottom
lattice_with_vacuum = Lattice.from_parameters(*lattice_params_vacuum, *angles_vacuum)

# Create vacuum layers on top and bottom of structure1 and structure2
structure1_with_vacuum = Structure(lattice_with_vacuum, species=structure1.species, coords=structure1.cart_coords)
structure2_with_vacuum = Structure(lattice_with_vacuum, species=structure2.species, coords=structure2.cart_coords)

# Calculate the translation vector for structure2_with_vacuum
translation = np.array([lattice_params_vacuum[0] * (3/4), 0, 0])

# Apply the translation to the coordinates of structure2_with_vacuum
for site in structure2_with_vacuum:
    site.translate_sites([0], translation)

# Append atoms from structure2_with_vacuum to structure1_with_vacuum
for site in structure2_with_vacuum:
    structure1_with_vacuum.append(site.specie, site.coords, coords_are_cartesian=True)

# Save the combined structure to a CIF file
structure1_with_vacuum.to(filename="solid_fluid_interface3.cif")
print("Combined structure has been made")

