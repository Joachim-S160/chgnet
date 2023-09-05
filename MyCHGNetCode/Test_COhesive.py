from ase import Atoms
from ase.io import read, write

# Load your crystal structure from a file or create it programmatically
# For this example, let's assume you have a CIF file named "my_structure.cif"
structure = read("chgnet/MyCHGNetCode/cif_files_halides/TiBr4_mp-569814_computed.cif")

# Identify the molecule you want to extract
# You may need to know the atom indices of the atoms in the molecule
# For example, if you want the first molecule in the structure, you can do:
molecule_indices = 0, 1, 2, 3  # Replace with the actual atom indices of the molecule

# Create a new ASE Atoms object containing the molecule
molecule = Atoms(symbols=structure.get_chemical_symbols()[*molecule_indices],
                 positions=structure.get_positions()[*molecule_indices])

# Now, you have your molecule as an ASE Atoms object
# You can access its properties, such as coordinates and chemical composition
print("Molecule:")
print("Chemical Symbols:", molecule.get_chemical_symbols())
print("Coordinates (Angstrom):")
print(molecule.get_positions())

# If you need to save the molecule to a file, you can do:
write("my_molecule.xyz", molecule)  # Save as an XYZ file, for example
