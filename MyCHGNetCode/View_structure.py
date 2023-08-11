"""

Desription:
Transform trajectory data to a small number of cif files, each cif file contains the structure of one of the frames of the trajectory.

"""

from plot import get_molecule_name
from ase.io.trajectory import Trajectory
import numpy as np

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
    for index, frame_index in enumerate(list_of_frames):
        tnof = len(list_of_frames)
        atoms = traj[frame_index]
        # ovito reads digit per digit, so make sure to pad with zeros
        power = np.floor(np.log10(tnof))
        if index < 10**power:
            frame_index_new = f"{'0' * int(power)}{index}"
            print(frame_index_new)
        else: frame_index_new = frame_index
        cif_filename = f'chgnet/MyCHGNetCode/cif_files_frames/WCl6_NVT/{get_molecule_name(trajfile)}_frame_{frame_index_new}.cif'
        atoms.write(cif_filename, format='cif')
    
def files_to_cif(files: list, number_of_frames: int=1):
    """
    input: 
        list of trajectory files: list of *.traj files
        number of frames: number of frames you want to convert to cif
    output: 
        cif files: *.cif
    """
    # Convert all trajectory files to CIF
    for file in files:
        traj = Trajectory(file)
        # tnof = total number of frames
        tnof = len(traj)
        print(f"Total number of frames in {file} is {tnof}")
        traj_to_cif(file, [int(tnof/number_of_frames * i) for i in range(number_of_frames)] + [tnof - 1])



# test
# files_to_cif(["chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj", "chgnet/MyCHGNetCode/mdNPT2_out_LiCl.traj", "chgnet/MyCHGNetCode/mdNPT2_out_TiBr4.traj", "chgnet/MyCHGNetCode/mdNPT2_out_TiI4.traj", "chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj", "chgnet/MyCHGNetCode/mdNPT3_out_Al.traj" ], 300)
files_to_cif(['chgnet/MyCHGNetCode/mdNVT_out_WCl6.traj'], 50)
print("Done!")
    