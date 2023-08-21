from pymatgen.core import Structure, Lattice
from pymatgen.io.vasp.outputs import Poscar
import numpy as np
import os

# # Load your structures from files or create them as needed
# structure1 = Structure.from_file("TiBr4.cif")
# structure2 = Structure.from_file("Tl2Cl3FAKE.cif")

# # Calculate the lattice parameters for the combined structure (same x lattice parameter, 4 times larger in y and z)
# lattice_params_vacuum = [structure1.lattice.abc[0], structure1.lattice.abc[1], structure1.lattice.abc[2]]
# print(f"Boxsize combined structure = {lattice_params_vacuum}")

# # Calculate the angles of the lattice
# angles_vacuum = structure1.lattice.angles

# # Create a new lattice with vacuum layers on top and bottom
# lattice_with_vacuum = Lattice.from_parameters(*lattice_params_vacuum, *angles_vacuum)

# # Create vacuum layers on top and bottom of structure1 and structure2
# structure1_with_vacuum = Structure(lattice_with_vacuum, species=structure1.species, coords=structure1.cart_coords)
# structure2_with_vacuum = Structure(lattice_with_vacuum, species=structure2.species, coords=structure2.cart_coords)

# # Calculate the translation vector for structure2_with_vacuum
# translation = np.array([lattice_params_vacuum[0] * (3/4), 0, 0])

# # Apply the translation to the coordinates of structure2_with_vacuum
# for site in structure2_with_vacuum:
#     site.translate_sites([0], translation)

# # Append atoms from structure2_with_vacuum to structure1_with_vacuum
# for site in structure2_with_vacuum:
#     structure1_with_vacuum.append(site.specie, site.coords, coords_are_cartesian=True)

# # Save the combined structure to a CIF file
# structure1_with_vacuum.to(filename="solid_fluid_interface3.cif")
# print("Combined structure has been made")

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

Biggest_elongated_box(structure = Structure.from_file("chgnet/MyCHGNetCode/cif_files_halides/LiCl_mp-22905_computed.cif"))

def test_Biggest_elongated_box():
    # Biggest_elongated_box(structure = Structure.from_file("chgnet/MyCHGNetCode/cif_files_halides/GeBr4_mp-567604_computed.cif"))
    

    folder_path = 'chgnet/MyCHGNetCode/cif_files_halides'
    file_extension = '.cif'  # Change this to the desired file extension

    for filename in os.listdir(folder_path):
        if filename.endswith(file_extension) and os.path.isfile(os.path.join(folder_path, filename)):
            print("File:", filename)
            Biggest_elongated_box(structure = Structure.from_file(os.path.join(folder_path, filename)))
    

# test_Biggest_elongated_box()