"""

Description:
2 phase coexistence model data processing:
detecting the melting point and error estimation

"""

import matplotlib.pyplot as plt
import numpy as np
from get_Data import *
from GetName import get_molecule_name

# use different colors for different files
colors = ['ro', 'bo', 'go', 'yo', 'mo', 'co', 'ko',
                     'rs', 'bs', 'gs', 'ys', 'ms', 'cs', 'ks',
                     'r^', 'b^', 'g^', 'y^', 'm^', 'c^', 'k^']

def time_temperature_plot(files: list, window_size: int=0.7):
    """
    Args: 
        list of trajectory files
        window_size: window size for time averaging
    returns: plot of time vs temperature    
    """
    for file in files:
        time = get_time(file)/2 # convert from 2 fs timestep to 1 fs timestep
        temperature = get_temperature(file)
        name = get_molecule_name(file)
        plt.plot(time, temperature, colors[files.index(file)], label= name)
        # Plot average temperature point with error bars on the plot
        average_temp = time_averaging(temperature, window_size)
        last_time = time[-1]
        y_err = convergence_error(temperature, window_size)
        plt.errorbar(last_time, average_temp, yerr=y_err, fmt='o', color='black', label='Melting point '+ name +': ' +
                 str(np.round(average_temp, 2)) + ' +/- ' + str(round(y_err, 2)) + ' K')  
        plt.xlabel('time (ps)')
        plt.ylabel('temperature (K)')
    plt.legend()
    plt.show()
time_temperature_plot(["chgnet/MyCHGNetCode/data_out_2PC/mdNVE_out_LiCl_cr_desta.traj"])

# Do the same for pressure:

def time_pressure_plot(files: list, window_size: int=0.7):
    """
    Args:   
        list of trajectory files
        window_size: window size for time averaging
    returns: plot of time vs pressure    
    """
    for file in files:
        time = get_time(file)/2 # convert from 2 fs timestep to 1 fs timestep
        pressure = get_pressure(file)
        name = get_molecule_name(file)
        plt.plot(time, pressure, colors[files.index(file)], label= name)
        # Plot average pressure point with error bars on the plot
        average_pressure = time_averaging(pressure, window_size)
        last_time = time[-1]
        y_err = convergence_error(pressure, window_size)
        plt.errorbar(last_time, average_pressure, yerr=y_err, fmt='o', color='black', label='Melting point '+ name +': ' +
                 str(np.round(average_pressure, 2)) + ' +/- ' + str(round(y_err, 2)) + ' GPa')  
        plt.xlabel('time (ps)')
        plt.ylabel('pressure (GPa)')
    plt.legend()
    plt.show()