"""

Desription:
Transform trajectory data (ASE) to a small number of cif files, each cif file contains the structure of one of the frames of the trajectory.

"""

from plot import get_molecule_name
from ase.io.trajectory import Trajectory
import numpy as np

def traj_to_cif(trajfile:str="chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj", stor_dir:str="chgnet/MyCHGNetCode/cif_files_frames/WCl6_NPT1", list_of_frames:list=[0]):
    """
    Converts selected frames from a trajectory file to CIF format.

    Args:
        trajfile (str): Path to the trajectory file (*.traj).
        stor_dir (str): Directory where the CIF files will be stored.
        list_of_frames (list): List of frame indices to convert to CIF.

    Returns:
        None. CIF files are saved to the specified directory.
    """
    # Load the trajectory file
    traj = Trajectory(trajfile)

    # Convert selected frames to CIF
    for index, frame_index in enumerate(list_of_frames):
        index +=1
        tnof = len(list_of_frames)
        atoms = traj[frame_index]
        # ovito reads digit per digit, so make sure to pad with zeros
        power = np.floor(np.log10(tnof))
        index_power = np.floor(np.log10(index))
        if index_power <= power:
            frame_index_new = f"{'0' * int(power - index_power)}{index}"
        else: 
            frame_index_new = frame_index
        cif_filename = f'{stor_dir}/{get_molecule_name(trajfile)}_frame_{frame_index_new}.cif'
        atoms.write(cif_filename, format='cif')
    
def files_to_cif(files: list, stor_dir: str, number_of_frames: int=1):
    """
    input: 
        list of trajectory files: list of *.traj files
        stor_dir: directory where the cif files are stored
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
        traj_to_cif(file, stor_dir, [int(tnof/number_of_frames * i) for i in range(number_of_frames)] + [tnof - 1])



# test
# files_to_cif(["chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj", "chgnet/MyCHGNetCode/mdNPT2_out_LiCl.traj", "chgnet/MyCHGNetCode/mdNPT2_out_TiBr4.traj", "chgnet/MyCHGNetCode/mdNPT2_out_TiI4.traj", "chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj", "chgnet/MyCHGNetCode/mdNPT3_out_Al.traj" ], 300)
# traj_to_cif(trajfile='chgnet/MyCHGNetCode/2PCTrajectories/mdNVT_out_Al.traj', stor_dir='chgnet/MyCHGNetCode/2PC_Cif_Files', list_of_frames=[-1])
files_to_cif(["chgnet/MyCHGNetCode/2PCTrajectories/mdNVT_out_Al.traj"], "chgnet/MyCHGNetCode/2PC_Cif_Files", 50)
print("Done!")
    