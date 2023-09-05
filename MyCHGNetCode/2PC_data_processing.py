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

def time_total_energy_plot(files:list) -> None:
    """
    Args: list of trajectory files
    returns: time vs total energy plot    
    """
    for file in files:
        time = get_time(file)/2 # convert from 2 fs timestep to 1 fs timestep
        Etot = get_total_energy(file)
        plt.plot(time, Etot, colors[files.index(file)], marker='o', linestyle='none', label= f"{get_molecule_name(file)} total energy")
    plt.xlabel('time (ps)')
    plt.ylabel('energy (eV)')
    plt.legend()
    plt.show()

# Average based on decimal cut-off
def discrete_averaging_system(x, y):
    """
    Args:
        x: x-axis values
        y: y-axis values
    Returns:
        unique x values with corresponding averaged y values for every unit of x
    """
    x_uniques, indices = np.unique(
        np.around(x, decimals=0), return_inverse=True)
    AV = np.zeros(len(x_uniques))
    for index_x_uniques, x_unique in enumerate(x_uniques):
        temporary = []
        for index_x, x_value in enumerate(x):
            if (x_value <= x_unique + 0.5) and (x_value >= x_unique - 0.5):
                temporary.append(y[index_x])
        AV[index_x_uniques] = np.average(temporary)
    return x_uniques, AV

# Average out quantity every 1 unit or more
def discrete_averaged_quantity(x, y, averagingfactor=1):
    """
    Args:
        x: x-axis values
        y: y-axis values
        averagingfactor: averaging factor for unit averaging
    """
    x /= averagingfactor
    x_, y_ = discrete_averaging_system(x, y)
    return x_*averagingfactor, y_

def time_averaged_Etot_plot(files: list, averagingfactor: int=1):
    """
    Args: 
        files: list of trajectory files
        averagingfactor: averaging factor for time averaging
    Returns:
        plot of time vs discrete averaged (total energy - time average)
    """
    
    for file in files:
        time = get_time(file)[:]/2 # convert from 2 fs timestep to 1 fs timestep
        number_of_atoms = get_number_of_atoms(file)[:]
        Etot = get_total_energy(file)[:]
        # Change to per atom
        Etot_per_atom = Etot/number_of_atoms
        plt.plot(time, Etot_per_atom, [color for color in reversed(colors)][files.index(file)], marker='o', linestyle='none', label= f"{get_molecule_name(file)} total energy per atom")
        # Do unit averagin
        discrete_averaged_time, discrete_averaged_Etot = discrete_averaged_quantity(time, Etot_per_atom, averagingfactor)
        plt.plot(discrete_averaged_time, discrete_averaged_Etot, colors[files.index(file)], marker='o', linestyle='none', label= f"{get_molecule_name(file)} total energy per atom, unit averaged")
    plt.xlabel('time (ps)')
    plt.ylabel('Total energy (eV)')
    plt.legend()
    plt.show()

def time_temperature_plot(files: list,extra_tags: list= [""], window_size: int=0.7 ):
    """
    Args: 
        files: list of trajectory files
        extra_tags: extra name tags to add
        window_size: window size for time averaging
    returns: plot of time vs temperature    
    """
    for file, extra_tag in zip(files, extra_tags):
        time = get_time(file)/2 # convert from 2 fs timestep to 1 fs timestep
        temperature = get_temperature(file)
        name = get_molecule_name(file, extra_tag=extra_tag)
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
# time_temperature_plot(["chgnet/MyCHGNetCode/data_out_2PC/mdNVE_out_LiCl_cr_desta.traj"])
# time_temperature_plot(["chgnet/MyCHGNetCode/data_out_2PC/mdNVE_out_LiCl_combined_relaxed_expanded.traj"])
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
        plt.errorbar(last_time, average_pressure, yerr=y_err, fmt='o', color='black', label='Pressure Melting Point '+ name +': ' +
                 str(np.round(average_pressure, 2)) + ' +/- ' + str(round(y_err, 2)) + ' GPa')  
        plt.xlabel('time (ps)')
        plt.ylabel('pressure (GPa)')
    plt.legend()
    plt.show()


# time_temperature_plot(["chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_LiBr_junction_relaxed.traj",
#                        "chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_LiBr_junction_relaxed_800K.traj",
#                        "chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_LiBr_junction_relaxed_1000K.traj"],["_2000K","_800K", "_1000K"])
time_temperature_plot(["chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_LiF_junction_relaxed.traj",
                       "chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_LiF_junction_relaxed_2250K.traj"],["_2000K","_2250K"])

# time_pressure_plot(["chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_Li2O_junction_relaxed.traj",
                    # "chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_Li2S_junction_relaxed.traj",
                    # "chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_Li2Se_junction_relaxed.traj"])
# time_averaged_Etot_plot(["chgnet/MyCHGNetCode/data_out_2PC/mdNVE_out_LiCl_combined_relaxed_expanded.traj"])
# time_temperature_plot(["chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_Li2O_junction_relaxed.traj",
                    # "chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_Li2S_junction_relaxed.traj",
                    # "chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_Li2Se_junction_relaxed.traj"])
# time_temperature_plot(["chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_LiF_junction_relaxed.traj"])
# time_averaged_Etot_plot(["chgnet/MyCHGNetCode/data_out_2PC_pymatgenjunction/mdNVE_out_LiF_junction_relaxed.traj"])
# time_pressure_plot(["chgnet/MyCHGNetCode/data_out_2PC/mdNVE_out_Al_cr_desta.traj"])
# time_total_energy_plot(["chgnet/MyCHGNetCode/data_out_2PC/mdNVE_out_Al_cr_desta.traj"])
# time_total_energy_plot(["chgnet/MyCHGNetCode/data_out_2PC/mdNVE_out_LiCl_cr_desta.traj"])
# time_averaged_Etot_plot(["chgnet/MyCHGNetCode/data_out_2PC/mdNVE_out_LiCl_cr_desta.traj"])
