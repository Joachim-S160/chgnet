"""

Description:
Plotting software for CHGNet data
Multiple files can be plotted on the same graph

"""

import matplotlib.pyplot as plt
import numpy as np
from get_Data import get_THtrp, get_tKEPE, normalize



def get_molecule_name(file:str):
    """
    Args: file
    returns: molecule name
    """
    
    input_string = file

    # Get the index of the last occurrence of '/'
    last_slash_index = input_string.rfind('/')

    # Get the index of the first occurrence of '.'
    dot_index = input_string.find('.')

    # Extract the substring between the last slash and the first dot (exclusive)
    desired_substring = input_string[last_slash_index + 11:dot_index]

    return desired_substring


def time_temperature_plot(files: list):
    """
    Args: list of trajectory files
    returns: plot of time vs temperature    
    """
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
    """
    Args: list of trajectory files
    returns: plot of time vs density    
    """
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
    
def temperature_enthalpy_plot(files: list, FP:bool=True):
    """
    Args: list of trajectory files
            FP: Fixed pressure; if True, pressure is fixed at 1 atm
    returns: plot of temperature vs enthalpy    
    """
    # use different colors for different files
    colors = ['ro', 'bo', 'go', 'yo', 'mo', 'co', 'ko', 'wo']
    THtr = get_THtrp(files, Fixed_pressure=FP)
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
    """
    Args: list of trajectory files
    returns: plot of temperature vs density    
    """
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
    """
    Args: list of trajectory files
    returns: time vs pressure plot
    """
    # use different colors for different files
    colors = ['ro', 'bo', 'go', 'yo', 'mo', 'co', 'ko', 'wo']
    THtr = get_THtrp(files)
    for index in range(len(files)):
        time = THtr[index, 2]
        pressure = THtr[index, 4]
        
        plt.plot(time, pressure/1000, colors[index], label= get_molecule_name(files[index]))
        plt.xlabel('time (ps)')
        plt.ylabel('pressure (atm)')
        # print(f"Linear regression for {files[index]} = {linear_regression(time, pressure)} in GPa/ps")
    plt.legend()
    plt.show()

def time_kinetic_energy_and_potential_energy_plot(files: list):
    """
    Args: list of trajectory files
    returns: time vs kinetic energy and potential energy plot    
    """
    # use different colors for different files
    colors = ['r', 'b', 'g', 'y', 'm', 'c', 'k', 'w']
    KEPEfiles = get_tKEPE(files)
    for index in range(len(files)):
        time = KEPEfiles[index, 0]
        kinetic_energy = KEPEfiles[index, 1]
        potential_energy = KEPEfiles[index, 2]
        
        plt.plot(time, normalize(kinetic_energy), colors[index], marker='o', linestyle='none', label= get_molecule_name(files[index]) + ' kinetic energy')
        plt.plot(time, normalize(potential_energy), colors[index], marker= 'x', linestyle='none', label= get_molecule_name(files[index]) + ' potential energy')
        plt.xlabel('time (ps)')
        plt.ylabel('energy (eV)')
        # print(f"Linear regression for {files[index]} = {linear_regression(time, kinetic_energy)} in eV/ps")
        # print(f"Linear regression for {files[index]} = {linear_regression(time, potential_energy)} in eV/ps")
    plt.legend()
    plt.show()

def linear_regression(x, y):
    """
    Args: x and y values
    returns: slope and intercept of linear regression
    """
    a, b = np.polyfit(x, y, deg=1)
    return a, b

def linear_fit(x, y):
    """
    Args: x and y values
    returns: slope and intercept of linear regression
    """
    A = np.vstack([x, np.ones(len(x))]).T
    a, b = np.linalg.lstsq(A, y, rcond=None)[0]
    return a, b



# density temperature discrete points diff md's, equilibrating each for 50 ps, WCl6, default taut, 100K steps around melting point, 6 points

# time_kinetic_energy_and_potential_energy_plot(['chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj','chgnet/MyCHGNetCode/mdNPT2_out_LiCl.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiBr4.traj','chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiI4.traj'])
# temperature_enthalpy_plot(['chgnet/MyCHGNetCode/mdNPT3_out_Al.traj'])


time_pressure_plot(['chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj'])
# temperature_enthalpy_plot(['chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj'])
# temperature_enthalpy_plot(['chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj','chgnet/MyCHGNetCode/mdNPT2_out_LiCl.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiBr4.traj','chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiI4.traj'])
# time_temperature_plot(['chgnet/MyCHGNetCode/mdNPT2_out_TiBr4.traj','chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiI4.traj'])
# time_temperature_plot(['chgnet/MyCHGNetCode/mdNPT2_out_LiCl.traj','chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj'])
# time_density_plot(['chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj','chgnet/MyCHGNetCode/mdNPT2_out_LiCl.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiBr4.traj','chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiI4.traj'])
# temperature_density_plot(['chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj','chgnet/MyCHGNetCode/mdNPT2_out_LiCl.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiBr4.traj','chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiI4.traj'])

def loopdyloop(listoffiles):
    """
    Args: list of files
    returns: plots    
    """
    for file in listoffiles:
        time_temperature_plot([file])
        temperature_density_plot([file])
        temperature_enthalpy_plot([file])
        
# loopdyloop(['chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj','chgnet/MyCHGNetCode/mdNPT2_out_LiCl.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiBr4.traj','chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiI4.traj'])