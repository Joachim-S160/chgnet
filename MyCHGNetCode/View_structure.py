"""

Desription:
Transform trajectory data to a small number of cif files, each cif file contains the structure of one of the frames of the trajectory.

"""

from plot import get_molecule_name
from ase.io.trajectory import Trajectory

def traj_to_cif(trajfile:str="chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj", list_of_frames:list=[0]):
    """
    input: 
        trajectory file: *.traj
        list of frames: list of the frames (indices) you want to convert to cif
    output: 
        cif file: *.cif
    """
    # Load the trajectory file
    traj = Trajectory(trajfile)

    # Convert selected frames to CIF
    for frame_index in list_of_frames:
        atoms = traj[frame_index]
        cif_filename = f'chgnet/MyCHGNetCode/cif_files_frames/frame_{frame_index}.cif'
        atoms.write(cif_filename, format='cif')
    
traj_to_cif()
    