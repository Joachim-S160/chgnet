import matplotlib.pyplot as plt
import numpy as np
from ase.io.trajectory import Trajectory


import numpy as np

def get_THtr(files: list):
    """
    input = list of trajectory files, 
    output = 3D Matrix: 
    1st dimension is the file, 
    2nd dimension is temperature, enthalpy, time, density 
    3rd dimension is the time step
    """

    # Create a 3D array to store the data
    num_files = len(files)
    num_timesteps = max(len(Trajectory(file)) for file in files)
    THtrfiles = np.zeros((num_files, 4, num_timesteps))

    for index, file in enumerate(files):
        traj = Trajectory(file)[:]
        num_timesteps_file = len(traj)
        Volume = np.array([atoms.get_volume() for atoms in traj])
        Epot = np.array([atoms.get_potential_energy() for atoms in traj])
        # NPT simulation so pressure is constant
        Enthalpy = Epot + Volume
        print(index, file)

        # Fill the 3D array with data
        THtrfiles[index, 0, :num_timesteps_file] = [atoms.get_temperature() for atoms in traj]
        THtrfiles[index, 1, :num_timesteps_file] = Enthalpy
        THtrfiles[index, 2, :num_timesteps_file] = np.arange(num_timesteps_file) * 0.2
        THtrfiles[index, 3, :num_timesteps_file] = 1 / Volume

    return THtrfiles

def get_molecule_name(file):
    input_string = file

    # Get the index of the last occurrence of '/'
    last_slash_index = input_string.rfind('/')

    # Get the index of the first occurrence of '.'
    dot_index = input_string.find('.')

    # Extract the substring between the last slash and the first dot (exclusive)
    desired_substring = input_string[last_slash_index + 12:dot_index]

    return desired_substring


def time_temperature_plot(files: list):
    # use different colors for different files
    colors = ['ro', 'bo', 'go', 'yo', 'mo', 'co', 'ko', 'wo']
    THtr = get_THtr(files)
    for index in range(len(files)):
        time = THtr[index, 2]
        temperature = THtr[index, 0]
        
        plt.plot(time, temperature, colors[index], label= get_molecule_name(files[index]))
        plt.xlabel('time (ps)')
        plt.ylabel('temperature (K)')
        print(f"Linear regression for {files[index]} = {linear_regression(time, temperature)}")
    plt.legend()
    plt.show()
    
def time_density_plot(files: list):
    # use different colors for different files
    colors = ['ro', 'bo', 'go', 'yo', 'mo', 'co', 'ko', 'wo']
    THtr = get_THtr(files)
    for index in range(len(files)):
        time = THtr[index, 2]
        density = THtr[index, 3]
        
        plt.plot(time, density, colors[index], label= get_molecule_name(files[index]))
        plt.xlabel('time (ps)')
        plt.ylabel('density (A^3)')
        print(f"Linear regression for {files[index]} = {linear_regression(time, density)}")
    plt.legend()
    plt.show()
    
def temperature_enthalpy_plot(files: list):
    # use different colors for different files
    colors = ['ro', 'bo', 'go', 'yo', 'mo', 'co', 'ko', 'wo']
    THtr = get_THtr(files)
    for index in range(len(files)):
        temperature = THtr[index, 0]
        enthalpy = THtr[index, 1]
        
        plt.plot(temperature, enthalpy, colors[index], label= get_molecule_name(files[index]))
        plt.xlabel('temperature (K)')
        plt.ylabel('enthalpy (eV)')
        print(f"Linear regression for {files[index]} = {linear_regression(temperature, enthalpy)}")
    plt.legend()
    plt.show()

def linear_regression(x, y):
    a, b = np.polyfit(x, y, deg=1)
    return a, b

def linear_fit(x, y):
    import numpy as np
    A = np.vstack([x, np.ones(len(x))]).T
    a, b = np.linalg.lstsq(A, y, rcond=None)[0]
    return a, b
  
# time_temperature_plot(['chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj','chgnet/MyCHGNetCode/mdNPT2_out_LiCl.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiBr4.traj','chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiI4.traj'])
time_density_plot(['chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj','chgnet/MyCHGNetCode/mdNPT2_out_LiCl.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiBr4.traj','chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiI4.traj'])