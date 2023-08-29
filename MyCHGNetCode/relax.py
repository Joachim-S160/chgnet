from __future__ import annotations
from chgnet.model import MolecularDynamics
from chgnet.model import StructOptimizer
from pymatgen.core import Structure
from chgnet.model import CHGNet


def relax(input_file: str = "POSCAR", name: str = "Al_combined", GPU: str = "cuda:0"):
    """
    Args: 
        input_file: input file of the structure
        name: name of the structure
        GPU: GPU to use
    Returns:
        None. Stores relaxed structure in {name}_relaxed.cif file
    """

    # load model
    chgnet = CHGNet.load()
    # load structure
    structure = Structure.from_file(input_file)
    print(f"Number of atoms: {structure.num_sites}")
    print('structure and model loaded')

    # Perform ionic relaxation (no Volume change)
    # = Relax the structure so that the atoms are moved to positions with lower potential energy but the cell size is NOT adjusted to the optimal size with no stresses.
    relaxer = StructOptimizer(use_device=GPU)
    relaxed_structure_dict: dict = relaxer.relax(
        structure, relax_cell=False)
    relaxed_structure: Structure = relaxed_structure_dict["final_structure"]
    print("Relaxation finished")
    relaxed_structure.to(filename=name + "_relaxed" + ".cif")
    print('relaxed structure saved')


relax("POSCAR", "Al_combined")
