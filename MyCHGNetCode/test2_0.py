# from pymatgen.core import Structure, Lattice
# from pymatgen.io.vasp.outputs import Poscar
# import numpy as np
# import os



# from pymatgen.core import Structure, Lattice
# from pymatgen.io.vasp.outputs import Poscar
# import numpy as np

# # Load the structures from files or create them as needed
# structure1 = Structure.from_file("chgnet/MyCHGNetCode/2PC_strcutures/LiCL/Fluid_LiCl.cif")
# structure2 = Structure.from_file("chgnet/MyCHGNetCode/2PC_strcutures/LiCL/Solid_LiCl.cif")

# borderlength = 2 #Angstrom
# # Calculate the lattice parameters for the vacuum layers (twice the original x lattice parameter)
# # assert structure2.lattice.abc == structure1.lattice.abc, "The two structures must have the same lattice parameters"
# lattice_params_vacuum = [2 * structure1.lattice.abc[0] + borderlength, structure1.lattice.abc[1], structure1.lattice.abc[2]] # +border in Angstrom for borders
# print(f"Boxsize combined structure = {lattice_params_vacuum}")

# # Calculate the angles of the lattice
# # assert structure1.lattice.angles == structure2.lattice.angles, "The two structures must have the same lattice angles"
# angles_vacuum = structure2.lattice.angles

# # Create a new lattice with vacuum layers 
# lattice_with_vacuum = Lattice.from_parameters(*lattice_params_vacuum, *angles_vacuum)

# # Create vacuum layers on top and bottom of structure1 and structure2
# border = np.array([[borderlength/2 , 0, 0] for _ in range(structure2.cart_coords.shape[0])])
# assert border.shape == structure2.cart_coords.shape, f"The border (shape = {border.shape}) vector must have the same shape as the coordinates of structure2 (shape = {structure2.cart_coords.shape}))"
# structure1_with_vacuum = Structure(lattice_with_vacuum, species=structure1.species, coords=structure1.cart_coords + border, coords_are_cartesian=True)
# translation = np.array([[lattice_params_vacuum[0]*(1/2), 0, 0] for _ in range(structure2.cart_coords.shape[0])])
# assert translation.shape == structure2.cart_coords.shape, "The translation vector must have the same shape as the coordinates of structure2"
# structure2_with_vacuum = Structure(lattice_with_vacuum, species=structure2.species, coords=structure2.cart_coords + border + translation, coords_are_cartesian=True)


# # Append atoms from structure2_with_vacuum to structure1_with_vacuum
# for site in structure2_with_vacuum:
#     structure1_with_vacuum.append(site.specie, site.frac_coords, coords_are_cartesian=False)

# # Save the combined structure to a CIF file
# structure1_with_vacuum.to(filename="solid_fluid_interface_Al.cif")
# print("combined structure has been made")

# # Melted to fast with the heating method: WCl6, TiI4, TiBr4,

# # def Biggest_elongated_box(structure=None, a=5, b=1, c=1, x_factor=5, Max_atoms=600):
# #     """
# #     Args:   
# #         structure: pymatgen structure 
# #         a,b,c: initial box dimensions or lattice parameters,
# #         x_factor: ratio between a and max(b,c) with a,b,c the lattice parameters
# #     Returns: 
# #         dimensions of the biggest possible supercell box  
    
# #     Note this function prefers the stretching of the box in the x direction
# #     Has extra functionality built in so that b and c lattice parameters are always bigger than 10 Angstrom
# #     """
# #     # Get number of atoms in unitcell

# #     noau = structure.num_sites
# #     A_lattice, B_lattice, C_lattice = structure.lattice.abc

# #     print(f"Number of atoms in unitcell: ", noau)

# #     noa = a * b * c * noau
# #     assert noa != 0, "No atoms in unitcell"
# #     assert noa < Max_atoms, f"Box is already bigger than {Max_atoms} atoms"
# #     # Update boxdimensions until noa is bigger than 500 or the boxdimensions are the same as before
# #     while noa < Max_atoms:
# #         x, y, z = a, b, c


# #         if a * (b + 1) * c * noau < Max_atoms and a >= x_factor * max(b + 1, c):
# #             b += 1

# #         if a * b * (c + 1) * noau < Max_atoms and a >= x_factor * max(b, c + 1):
# #             c += 1
# #         if (a + 1) * b * c * noau < Max_atoms:
# #             a += 1
# #         if [x, y, z] == [a, b, c]:
# #             break
# #         noa = a * b * c * noau
# #     print("Box dimensions: ", a, b, c)
# #     print("Total number of atoms in the box: ", a*b*c*noau)
# #     print(f"Boxsize supercell structure = {a * A_lattice, b * B_lattice, c * C_lattice} [Angstrom]")
# #     assert b * B_lattice >= 10, f"b lattice parameter is smaller than 10 Angstrom, b lattice parameter is {b * B_lattice}"
# #     assert c * C_lattice >= 10, f"c lattice parameter is smaller than 10 Angstrom, c lattice parameter is {c * C_lattice}"
# #     assert a *A_lattice * 5 >= b * B_lattice, f"Box is not elongated enough, a lattice parameter is {a * A_lattice} and b lattice parameter is {b * B_lattice}"
# #     assert a *A_lattice * 5 >= c * C_lattice, f"Box is not elongated enough, a lattice parameter is {a * A_lattice} and c lattice parameter is {c * C_lattice}"

# #     return [a, b, c]

# # Biggest_elongated_box(structure = Structure.from_file("chgnet/MyCHGNetCode/cif_files_halides/TiI4_mp-541013_computed.cif"))

# # def test_Biggest_elongated_box():
# #     # Biggest_elongated_box(structure = Structure.from_file("chgnet/MyCHGNetCode/cif_files_halides/GeBr4_mp-567604_computed.cif"))
    

# #     folder_path = 'chgnet/MyCHGNetCode/cif_files_halides'
# #     file_extension = '.cif'  # Change this to the desired file extension

# #     for filename in os.listdir(folder_path):
# #         if filename.endswith(file_extension) and os.path.isfile(os.path.join(folder_path, filename)):
# #             print("File:", filename)
# #             Biggest_elongated_box(structure = Structure.from_file(os.path.join(folder_path, filename)))
    

# # test_Biggest_elongated_box()

import math

def triclinic_to_cartesian(x_frac, y_frac, z_frac, a, b, c, alpha, beta, gamma):
    # Convert angles to radians
    alpha = math.radians(alpha)
    beta = math.radians(beta)
    gamma = math.radians(gamma)
    
    # Calculate unit cell vectors in Cartesian coordinates
    a1 = [a, 0, 0]
    a2 = [b * math.cos(gamma), b * math.sin(gamma), 0]
    a3 = [c * math.cos(beta), c * (math.cos(alpha) - math.cos(beta) * math.cos(gamma)) / math.sin(gamma), c * math.sqrt(1 - math.cos(beta)**2 - math.cos(alpha)**2 - math.cos(gamma)**2 + 2 * math.cos(alpha) * math.cos(beta) * math.cos(gamma)) / math.sin(gamma)]
    
    # Convert fractional coordinates to Cartesian
    x_cart = x_frac * a1[0] + y_frac * a2[0] + z_frac * a3[0]
    y_cart = x_frac * a1[1] + y_frac * a2[1] + z_frac * a3[1]
    z_cart = x_frac * a1[2] + y_frac * a2[2] + z_frac * a3[2]
    
    return x_cart, y_cart, z_cart

# Example usage
x_frac = 0
y_frac =1
z_frac = 0
a = 5.0
b = 4.0
c = 6.0
alpha = 60.0
beta = 60.0
gamma = 60.0

x_cart, y_cart, z_cart = triclinic_to_cartesian(x_frac, y_frac, z_frac, a, b, c, alpha, beta, gamma)
print(f"Cartesian coordinates: ({x_cart}, {y_cart}, {z_cart})")
