"""

The purpose of this python script is to calculate the cohesive energy of a structure using pymatgen

"""

from __future__ import annotations
from chgnet.model import StructOptimizer
from pymatgen.core import Structure, Lattice
from chgnet.model import CHGNet
import numpy as np
import matplotlib.pyplot as plt

# load model
chgnet = CHGNet.load()


def ionic_relax(struct: Structure) -> Structure:
    """
    Args: 
        struct: structure to be relaxed
    Returns:
        ionic relaxed structure
    """
    # Perform ionic relaxation (no Volume change)
    # = Relax the structure so that the atoms are moved to positions with lower potential energy but the cell size is NOT adjusted to the optimal size with no stresses.
    relaxer = StructOptimizer()
    relaxed_structure_dict: dict = relaxer.relax(
        struct, relax_cell=False)
    relaxed_struct: Structure = relaxed_structure_dict["final_structure"]
    print("Relaxation finished")
    return relaxed_struct


def get_chgnet_predicted_energy(struct: Structure) -> float:
    """
    Args:
        struct: structure
    Returns:
        predicted (by chgnet) energy of the structure
    """

    prediction = chgnet.predict_structure(struct)
    key, _unit = "energy", "eV/atom"
    return prediction[key[0]]


def Cohesive_energy(unit_cell: str = "LiCl.cif", single_molecule: str = "LiCl_molecule.cif", base_name: str = "LiCl"):
    """
    Args:
        unit_cell (str): path to the unit cell cif file
        single_molecule (str): path to the single molecule cif file
        base_name (str): base name for the output files
    Returns:
        Cohesive energy of the structure
    """

    # load unit cell
    unit_cell_structure = Structure.from_file(unit_cell)

    # Get single molecule
    single_molecule_structure = Structure.from_file(single_molecule)

    # Create lots of empty space with molecule in the middle
    lattice_params_vacuum = [30, 30, 30]  # in Angstrom
    angles_vacuum = single_molecule_structure.lattice.angles
    # Create a new lattice with vacuum layers
    lattice_with_vacuum = Lattice.from_parameters(
        *lattice_params_vacuum, *angles_vacuum)
    Spaceous_molecule = Structure(lattice_with_vacuum, species=single_molecule_structure.species,
                                  coords=single_molecule_structure.cart_coords, coords_are_cartesian=True)

    Spaceous_molecule.to(
        filename="chgnet/MyCHGNetCode/Cohesive_cifs/" + base_name + "Spaceous_molecule.cif")

    relaxed_molecule = ionic_relax(Spaceous_molecule)

    relaxed_molecule.to(
        filename="chgnet/MyCHGNetCode/Cohesive_cifs/" + base_name + "relaxed_molecule.cif")

    relaxed_unit_cell = ionic_relax(unit_cell_structure)

    E_unit_cell = get_chgnet_predicted_energy(relaxed_unit_cell)
    print(f"Predicted energy of {base_name} unit cell: {E_unit_cell}")

    E_molecule = get_chgnet_predicted_energy(relaxed_molecule)
    print(f"Predicted energy of {base_name} molecule: {E_molecule}")

    E_cohesive = E_unit_cell - E_molecule
    print(
        f"Normalized cohesive energy of {base_name}: {E_cohesive/E_molecule}")
    return E_cohesive/E_molecule


# Cohesive_energy("chgnet/MyCHGNetCode/cif_files_halides/TiBr4_mp-569814_computed.cif",
#                 "chgnet/MyCHGNetCode/Cohesive_cifs/TiBr4_molecule.cif", "TiBr4")
# Cohesive_energy("chgnet/MyCHGNetCode/Cohesive_cifs/WCl6_mp-571518_computed.cif",
#                 "chgnet/MyCHGNetCode/Cohesive_cifs/WCl6_molecule.cif", "WCl6")

def plot_cohesive_energy(unitcells: list = ["chgnet/MyCHGNetCode/cif_files_halides/TiBr4_mp-569814_computed.cif"], molecules: list = ["chgnet/MyCHGNetCode/Cohesive_cifs/TiBr4_molecule.cif"], labels: list = ["TiBr4"]):
    """
    Args:
        unitcells: list of unitcell cif files of material
        molecules: list of single molecule cif files of material
        labels: list of labels for the materials
    Returns:
        plot of cohesive energy vs material
    """
    cohesive_energies = []
    for unitcell, molecule in zip(unitcells, molecules):
        cohesive_energies.append(
            Cohesive_energy(unit_cell=unitcell, single_molecule=molecule))
    ys = cohesive_energies
    xs = labels
    for x,y in zip(xs,ys):

        label = "{:.6f}".format(y)

        plt.annotate(label, # this is the text
                    (x,y), # these are the coordinates to position the label
                    textcoords="offset points", # how to position the text
                    xytext=(0,10), # distance from text to points (x,y)
                    ha='center') # horizontal alignment can be left, right or center
    
    plt.plot(labels, cohesive_energies, marker='o', linestyle='none')
    plt.xlabel('Material')
    plt.ylabel('Cohesive energy (eV)')
    plt.show()
    
plot_cohesive_energy(["chgnet/MyCHGNetCode/cif_files_halides/TiBr4_mp-569814_computed.cif",
                      "chgnet/MyCHGNetCode/Cohesive_cifs/WCl6_mp-571518_computed.cif"],[
                          "chgnet/MyCHGNetCode/Cohesive_cifs/TiBr4_molecule.cif",
                          "chgnet/MyCHGNetCode/Cohesive_cifs/WCl6_molecule.cif"],[
                              "TiBr4","WCl6"
                          ])