import matplotlib.pyplot as plt
import numpy as np
from ase.io.trajectory import Trajectory


# normalize function
def normalize(x):
    return (x-np.min(x))/(np.max(x)-np.min(x))

def get_THtrp(files: list):
    """
    input = list of trajectory files, 
    output = 3D Matrix: 
    1st dimension is the file, 
    2nd dimension is temperature, enthalpy, time, density, pressure 
    3rd dimension is the time step
    """

    # Create a 3D array to store the data
    num_files = len(files)
    num_timesteps = max(len(Trajectory(file)) for file in files)
    THtrfiles = np.zeros((num_files, 5, num_timesteps))

    for index, file in enumerate(files):
        traj = Trajectory(file)[:]
        num_timesteps_file = len(traj)
        Volume = np.array([atoms.get_volume() for atoms in traj])
        Etot = np.array([atoms.get_total_energy() for atoms in traj])
        # stress tensors = sigma_ij at each frame
        # sigma_ij = sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_xz, sigma_yz
        # pressure = - 1/3 * (sigma_xx + sigma_yy + sigma_zz)
        stresstensors = [atoms.get_stress() for atoms in traj]
        pressure = np.array([-1/3 * (sigma[0] + sigma[1] + sigma[2]) for sigma in stresstensors]) 
        Enthalpy = normalize(Etot + pressure * Volume)
        mass = traj[0].get_masses().sum()
        print(index, file)

        # Fill the 3D array with data
        THtrfiles[index, 0, :num_timesteps_file] = [atoms.get_temperature() for atoms in traj]
        THtrfiles[index, 1, :num_timesteps_file] = Enthalpy
        THtrfiles[index, 2, :num_timesteps_file] = np.arange(num_timesteps_file) * 0.2
        # Convert density from amu/A^3 to g/cm^3 by multiplying by 1.66054
        THtrfiles[index, 3, :num_timesteps_file] = 1.66054 * mass / Volume
        THtrfiles[index, 4, :num_timesteps_file] = pressure

    return THtrfiles

def get_t_KE_and_PE(files: list):
    """
    input = list of trajectory files, 
    output = 3D Matrix: 
    1st dimension is the file, 
    2nd dimension is time, kinetic energy, potential energy
    3rd dimension is the time step
    """

    # Create a 3D array to store the data
    num_files = len(files)
    num_timesteps = max(len(Trajectory(file)) for file in files)
    KEPEfiles = np.zeros((num_files, 3, num_timesteps))

    for index, file in enumerate(files):
        traj = Trajectory(file)[:]
        num_timesteps_file = len(traj)
        Epot = np.array([atoms.get_potential_energy() for atoms in traj])
        Ekin = np.array([atoms.get_kinetic_energy() for atoms in traj])
        print(index, file)

        # Fill the 3D array with data
        KEPEfiles[index, 0, :num_timesteps_file] = np.arange(num_timesteps_file) * 0.2
        KEPEfiles[index, 1, :num_timesteps_file] = Ekin
        KEPEfiles[index, 2, :num_timesteps_file] = Epot

    return KEPEfiles

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
    THtr = get_THtrp(files)
    for index in range(len(files)):
        time = THtr[index, 2]
        temperature = THtr[index, 0]
        name = get_molecule_name(files[index])
        plt.plot(time, temperature, colors[index], label= name)
        plt.xlabel('time (ps)')
        plt.ylabel('temperature (K)')
        print(f"Linear regression for {name} = {linear_regression(time, temperature)} in K/ps")
    plt.legend()
    plt.show()
    
def time_density_plot(files: list):
    # use different colors for different files
    colors = ['ro', 'bo', 'go', 'yo', 'mo', 'co', 'ko', 'wo']
    THtr = get_THtrp(files)
    for index in range(len(files)):
        time = THtr[index, 2]
        density = THtr[index, 3]
        
        plt.plot(time, density, colors[index], label= get_molecule_name(files[index]))
        plt.xlabel('time (ps)')
        plt.ylabel('density (gm/cc)')
        # print(f"Linear regression for {files[index]} = {linear_regression(time, density)} in gm/cc/ps")
    plt.legend()
    plt.show()
    
def temperature_enthalpy_plot(files: list):
    # use different colors for different files
    colors = ['ro', 'bo', 'go', 'yo', 'mo', 'co', 'ko', 'wo']
    THtr = get_THtrp(files)
    for index in range(len(files)):
        temperature = THtr[index, 0]
        enthalpy = THtr[index, 1]
        
        plt.plot(temperature, enthalpy, colors[index], label= get_molecule_name(files[index]))
        plt.xlabel('temperature (K)')
        plt.ylabel('enthalpy (eV)')
        # print(f"Linear regression for {files[index]} = {linear_regression(temperature, enthalpy)} in eV/K")
    plt.legend()
    plt.show()

def temperature_density_plot(files: list):
    # use different colors for different files
    colors = ['ro', 'bo', 'go', 'yo', 'mo', 'co', 'ko', 'wo']
    THtr = get_THtrp(files)
    for index in range(len(files)):
        temperature = THtr[index, 0]
        density = THtr[index, 3]
        
        plt.plot(temperature, density, colors[index], label= get_molecule_name(files[index]))
        plt.xlabel('temperature (K)')
        plt.ylabel('density (gm/cc)')
        # print(f"Linear regression for {files[index]} = {linear_regression(temperature, density)} in gm/cc/K")
    plt.legend()
    plt.show()

def time_pressure_plot(files: list):
    # use different colors for different files
    colors = ['ro', 'bo', 'go', 'yo', 'mo', 'co', 'ko', 'wo']
    THtr = get_THtrp(files)
    for index in range(len(files)):
        time = THtr[index, 2]
        pressure = THtr[index, 4]
        
        plt.plot(time, pressure, colors[index], label= get_molecule_name(files[index]))
        plt.xlabel('time (ps)')
        plt.ylabel('pressure (GPa)')
        # print(f"Linear regression for {files[index]} = {linear_regression(time, pressure)} in GPa/ps")
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



# density temperature discrete points diff md's, equilibrating each for 50 ps, WCl6, default taut, 100K steps around melting point, 6 points
# Ask Bowen if I can visualize the structure's movement in ASE or ovito
  
time_pressure_plot(['chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj','chgnet/MyCHGNetCode/mdNPT2_out_LiCl.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiBr4.traj','chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiI4.traj'])
# temperature_enthalpy_plot(['chgnet/MyCHGNetCode/mdNPT2_out_LiCl.traj','chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj'])
# temperature_enthalpy_plot(['chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj','chgnet/MyCHGNetCode/mdNPT2_out_LiCl.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiBr4.traj','chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiI4.traj'])
# time_temperature_plot(['chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj','chgnet/MyCHGNetCode/mdNPT2_out_LiCl.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiBr4.traj','chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiI4.traj'])
# time_density_plot(['chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj','chgnet/MyCHGNetCode/mdNPT2_out_LiCl.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiBr4.traj','chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiI4.traj'])
# temperature_density_plot(['chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj','chgnet/MyCHGNetCode/mdNPT2_out_LiCl.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiBr4.traj','chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiI4.traj'])