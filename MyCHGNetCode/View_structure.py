"""

Description:
Transform trajectory data (ASE) to a small number of cif files, each cif file contains the structure of one of the frames of the trajectory.

"""

from GetName import get_molecule_name
from ase.io.trajectory import Trajectory
import numpy as np
import os

def create_directory(directory_path):
    """
    Args:
        directory_path: path to the directory to be created
    Returns:
        None. Creates directory if it does not exist yet    
    """
    if not os.path.exists(directory_path):
        os.makedirs(directory_path)
        print(f"Directory '{directory_path}' created.")
    else:
        print(f"Directory '{directory_path}' already exists.")

def traj_to_cif(trajfile:str="chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj", stor_dir:str="chgnet/MyCHGNetCode/cif_files_frames/WCl6_NPT1", list_of_frames:list=[0]):
    """
    Args:
        trajfile (str): Path to the trajectory file (*.traj).
        stor_dir (str): Directory where the CIF files will be stored.
        list_of_frames (list): List of frame indices to convert to CIF.

    Returns:
        None. CIF files are saved to the specified directory.
        
    Note:   Converts selected frames from a trajectory file to CIF format.
    """
    # Make directory
    create_directory(stor_dir)
    
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
    
def files_to_cif(files: list, stor_dir: list, number_of_frames: int=1):
    """
    input: 
        list of trajectory files: list of *.traj files
        stor_dir: list of directories where the cif files are stored
        number of frames: number of frames you want to convert to cif
    output: 
        cif files: *.cif
    Note: for ovito GIF or MP4 purposes
            Make sure to have the same number of directories as trajectory files and have them in the same order
    """
    assert len(files) == len(stor_dir), "Number of trajectory files and number of directories do not match!"
    # Convert all trajectory files to CIF
    for file, dir in zip(files, stor_dir):
        traj = Trajectory(file)
        # tnof = total number of frames
        tnof = len(traj)
        print(f"Total number of frames in {file} is {tnof}")
        traj_to_cif(file, dir, [int(tnof/number_of_frames * i) for i in range(number_of_frames)] + [tnof - 1])


# files_to_cif(["chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_LiCl_junction_relaxed.traj",
#               "chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_Li2CO3_junction_relaxed.traj",
#               "chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_Li2O_junction_relaxed.traj",
#               "chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_Li2S_junction_relaxed.traj",
#               "chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_Li2Se_junction_relaxed.traj",
#               "chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_Li3N_junction_relaxed.traj",
#               "chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_LiBr_junction_relaxed.traj",
#               "chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_LiF_junction_relaxed.traj"
#               ], 
#              ["chgnet/MyCHGNetCode/2PC_Cif_Files_frames_PMJ/LiCl_PMJ",
#               "chgnet/MyCHGNetCode/2PC_Cif_Files_frames_PMJ/Li2CO3_PMJ",
#               "chgnet/MyCHGNetCode/2PC_Cif_Files_frames_PMJ/Li2O_PMJ",
#               "chgnet/MyCHGNetCode/2PC_Cif_Files_frames_PMJ/Li2S_PMJ",
#               "chgnet/MyCHGNetCode/2PC_Cif_Files_frames_PMJ/Li2Se_PMJ",
#                 "chgnet/MyCHGNetCode/2PC_Cif_Files_frames_PMJ/Li3N_PMJ",
#                 "chgnet/MyCHGNetCode/2PC_Cif_Files_frames_PMJ/LiBr_PMJ",
#                 "chgnet/MyCHGNetCode/2PC_Cif_Files_frames_PMJ/LiF_PMJ"
#               ], 100)

# files_to_cif([
#                 "chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_Li3N_junction_relaxed_1000K.traj",
#                 "chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_LiBr_junction_relaxed_1000K.traj",
#                 "chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_Li2CO3_junction_relaxed_1000K.traj"],
#              [
#                 "chgnet/MyCHGNetCode/2PC_Cif_Files_frames_PMJ/Li3N_PMJ_1000K",
#                 "chgnet/MyCHGNetCode/2PC_Cif_Files_frames_PMJ/LiBr_PMJ_1000K",
#                 "chgnet/MyCHGNetCode/2PC_Cif_Files_frames_PMJ/Li2CO3_PMJ_1000K"], 100)

files_to_cif(["chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_LiF_junction_relaxed_2250K.traj",
              "chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_LiCl_junction_relaxed_2600K.traj",
              "chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_LiBr_junction_relaxed_800K.traj",
              "chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_Li3N_junction_relaxed_600K.traj",
              "chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_Li2Se_junction_relaxed_1000K.traj",
              "chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_Li2S_junction_relaxed_1000K.traj",
              "chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_Li2O_junction_relaxed_2600K.traj"
              ],
             ["chgnet/MyCHGNetCode/2PC_Cif_Files_frames_PMJ/LiF_2250K",
              "chgnet/MyCHGNetCode/2PC_Cif_Files_frames_PMJ/LiCl_2600K",
              "chgnet/MyCHGNetCode/2PC_Cif_Files_frames_PMJ/LiBr_800K",
              "chgnet/MyCHGNetCode/2PC_Cif_Files_frames_PMJ/Li3N_600K",
              "chgnet/MyCHGNetCode/2PC_Cif_Files_frames_PMJ/Li2Se_1000K",
              "chgnet/MyCHGNetCode/2PC_Cif_Files_frames_PMJ/Li2S_1000K",
              "chgnet/MyCHGNetCode/2PC_Cif_Files_frames_PMJ/Li2O_2600K"
              ], 200)

print("Done!")
    