"""

Description:
Plotting software for CHGNet data
Multiple files can be plotted on the same graph

"""

import matplotlib.pyplot as plt
import numpy as np
from get_Data import *
from GetName import get_molecule_name

# use different colors for different files
colors = ['ro', 'bo', 'go', 'yo', 'mo', 'co', 'ko',
                     'rs', 'bs', 'gs', 'ys', 'ms', 'cs', 'ks',
                     'r^', 'b^', 'g^', 'y^', 'm^', 'c^', 'k^']

def time_temperature_plot(files: list):
    """
    Args: list of trajectory files
    returns: plot of time vs temperature    
    """
    for file in files:
        time = get_time(file)
        temperature = get_temperature(file)
        name = get_molecule_name(file)
        plt.plot(time, temperature, colors[files.index(file)], label= name)
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
    for file in files:
        time = get_time(file)
        density = get_density(file)
        name = get_molecule_name(file)
        plt.plot(time, density, colors[files.index(file)], label= name)
        plt.xlabel('time (ps)')
        plt.ylabel('density (gm/cc)')
    plt.legend()
    plt.show()
        
def temperature_enthalpy_plot(files: list):
    """
    Args: list of trajectory files
    returns: plot of temperature vs enthalpy    
    """
    for file in files:
        temperature = get_temperature(file)
        enthalpy = get_enthalpy(file)
        name = get_molecule_name(file)
        plt.plot(temperature, enthalpy, colors[files.index(file)], label= name)
        plt.xlabel('temperature (K)')
        plt.ylabel('enthalpy (eV)')
        print(f"Linear regression for {name} = {linear_regression(temperature, enthalpy)} in eV/K")
    plt.legend()
    plt.show()
    
def temperature_density_plot(files: list):
    """
    Args: list of trajectory files
    returns: plot of temperature vs density    
    """
    for file in files:
        temperature = get_temperature(file)
        density = get_density(file)
        name = get_molecule_name(file)
        plt.plot(temperature, density, colors[files.index(file)], label= name)
        plt.xlabel('temperature (K)')
        plt.ylabel('density (gm/cc)')
    plt.legend()
    plt.show()

def time_pressure_plot(files: list):
    """
    Args: list of trajectory files
    returns: time vs pressure plot
    """
    for file in files:
        time = get_time(file)
        pressure = get_pressure(file)
        name = get_molecule_name(file)
        plt.plot(time, pressure, colors[files.index(file)], label= name)
        plt.xlabel('time (ps)')
        plt.ylabel('pressure (GPa)')
    plt.legend()
    plt.show()

def time_kinetic_energy_and_potential_energy_plot(files: list):
    """
    Args: list of trajectory files
    returns: time vs kinetic energy and potential energy plot    
    """
    for file in files:
        time = get_time(file)
        # Ekin = get_kinetic_energy(file)
        # plt.plot(time, Ekin, colors[files.index(file)], marker='o', linestyle='none', label= f"{get_molecule_name(file)} kinetic energy")
        Epot = get_potential_energy(file)
        plt.plot(time, Epot, colors[files.index(file)], marker='x', linestyle='none', label= f"{get_molecule_name(file)} potential energy")
    plt.xlabel('time (ps)')
    plt.ylabel('energy (eV)')
    plt.legend()
    plt.show()

def time_total_energy_plot(files:list) -> None:
    """
    Args: list of trajectory files
    returns: time vs total energy plot    
    """
    for file in files:
        time = get_time(file)
        Etot = get_total_energy(file)
        plt.plot(time, Etot, colors[files.index(file)], marker='o', linestyle='none', label= f"{get_molecule_name(file)} total energy")
    plt.xlabel('time (ps)')
    plt.ylabel('energy (eV)')
    plt.legend()
    plt.show()

# For discrete temperature method only
def discrete_temperature_density_plot() -> None:
    """
    Args: None
    returns: None, plot of temperature vs density discrete points
    """
    for index, temp in enumerate(range(100,500,100)):
        temperature = get_temperature(f"chgnet/MyCHGNetCode/data_out_discrete_Tr/mdNPT_out_WCl6_{temp}.traj")
        density = get_density(f"chgnet/MyCHGNetCode/data_out_discrete_Tr/mdNPT_out_WCl6_{temp}.traj")
        plt.plot(time_averaging(temperature, window_factor = 0.95), time_averaging(density), colors[index], label= f"{temp} K")
    plt.xlabel('temperature (K)')
    plt.ylabel('density (gm/cc)')
    plt.legend()
    plt.show()
    
def density_time_plot_DISCRETE() -> None:
    """
    Args: None
    returns: None, plot of density vs time discrete points
    """
    for index, temp in enumerate(range(100,1000,100)):
        plt.plot(get_time(f"chgnet/MyCHGNetCode/data_out_discrete_Tr/mdNPT_out_WCl6_{temp}.traj"), get_density(f"chgnet/MyCHGNetCode/data_out_discrete_Tr/mdNPT_out_WCl6_{temp}.traj"), colors[index], label= f"{temp} K")
        plt.xlabel('time (ps)')
        plt.ylabel('density (gm/cc)')
        plt.legend()
        plt.show()
    

def linear_regression(x, y):
    """
    Args: x and y values
    returns: slope and intercept of linear regression
    """
    a, b = np.polyfit(x, y, deg=1)
    return a, b

time_total_energy_plot(["chgnet/MyCHGNetCode/data_out_discrete_Tr/mdNPT_out_Al_100.traj"])
# time_density_plot(["chgnet/MyCHGNetCode/data_out_heating/mdNPT2_out_LiCl.traj"])
# temperature_density_plot(["chgnet/MyCHGNetCode/data_out_heating/mdNPT2_out_LiCl.traj"])
# temperature_density_plot(["chgnet/MyCHGNetCode/data_out_2PC/mdNVE_out_LiCl_cr_desta.traj"])
# time_pressure_plot(["chgnet/MyCHGNetCode/data_out_discrete_Tr/mdNPT_out_Al_100.traj"])
# discrete_temperature_density_plot()
# density_time_plot_DISCRETE()
# density temperature discrete points diff md's, equilibrating each for 50 ps, WCl6, default taut, 100K steps around melting point, 6 points
# time_total_energy_plot(["chgnet/MyCHGNetCode/data_out_heating/mdNPT2_out_HfF4.traj"])
# time_kinetic_energy_and_potential_energy_plot(['chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj','chgnet/MyCHGNetCode/mdNPT2_out_LiCl.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiBr4.traj','chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiI4.traj'])
# temperature_enthalpy_plot(['chgnet/MyCHGNetCode/mdNPT3_out_Al.traj'])
# time_density_plot(["chgnet/MyCHGNetCode/data_out_discrete_Tr/mdNPT_out_WCl6_100.traj"])
# time_density_plot(["chgnet/MyCHGNetCode/data_out_discrete_Tr/mdNPT_out_WCl6_200.traj"])
# time_density_plot(["chgnet/MyCHGNetCode/data_out_discrete_Tr/mdNPT_out_WCl6_300.traj"])
# time_total_energy_plot(["chgnet/MyCHGNetCode/data_out_2PC/mdNVE_out_Al_cr_desta.traj", "chgnet/MyCHGNetCode/data_out_2PC/mdNVE_out_LiCl_cr_desta.traj"])
# time_pressure_plot(["chgnet/MyCHGNetCode/data_out_2PC/mdNVE_out_Al_cr_desta.traj", "chgnet/MyCHGNetCode/data_out_2PC/mdNVE_out_LiCl_cr_desta.traj"])
# time_density_plot(["chgnet/MyCHGNetCode/data_out_discrete_Tr/mdNPT_out_WCl6_100.traj"])
# temperature_density_plot(['chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj'])
# temperature_enthalpy_plot(['chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj'])
# temperature_enthalpy_plot(['chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj','chgnet/MyCHGNetCode/mdNPT2_out_LiCl.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiBr4.traj','chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiI4.traj'])
# time_temperature_plot(['chgnet/MyCHGNetCode/mdNPT2_out_TiBr4.traj','chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiI4.traj'])
# time_temperature_plot(['chgnet/MyCHGNetCode/mdNPT2_out_LiCl.traj','chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj'])
# time_density_plot(['chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj','chgnet/MyCHGNetCode/mdNPT2_out_LiCl.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiBr4.traj','chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiI4.traj'])
# temperature_density_plot(['chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj','chgnet/MyCHGNetCode/mdNPT2_out_LiCl.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiBr4.traj','chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiI4.traj'])
# time_temperature_plot(['chgnet/MyCHGNetCode/data_out_heating/mdNPT2_out_HfF4.traj','chgnet/MyCHGNetCode/data_out_heating/mdNPT2_out_TiBr4.traj', 'chgnet/MyCHGNetCode/data_out_heating/mdNPT2_out_TiI4.traj'])
# time_pressure_plot(["chgnet/MyCHGNetCode/data_out_heating/mdNPT2_out_HfF4.traj"])
# time_kinetic_energy_and_potential_energy_plot(["chgnet/MyCHGNetCode/data_out_heating/mdNPT2_out_HfF4.traj"])
# time_density_plot(["chgnet/MyCHGNetCode/data_out_heating/mdNPT1_out_TiI4.traj"])
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