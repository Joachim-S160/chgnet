from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
import heapq
from scipy.interpolate import make_interp_spline
from ase.io.trajectory import Trajectory

traj = Trajectory("chgnet/MyCHGNetCode/mdNPT3_out.traj")[:]

Volume = np.array([atoms.get_volume() for atoms in traj])
Epot = np.array([atoms.get_potential_energy() for atoms in traj])
Temperature = np.array([atoms.get_temperature() for atoms in traj])

# NPT simulation so pressure is constant
Enthalpy = Epot + Volume

# test out raw data
plt.plot(Temperature, Enthalpy, 'ro')
plt.show()

def plot_density2(txtfile):
    #  10 different colors
    data = np.loadtxt(txtfile, skiprows=1)

    time = data[:, 0]
    density = data[:, 1]

    plt.plot(time, density, color='g', label=f'density at 1200 K')
    plt.legend(loc='best')
    plt.xlabel('time [femtoseconds = 10^-15 s]]')
    plt.ylabel('density [g/cm^3]')
    plt.show()

################################### HEATING ######################################################################################

# normalize function
def normalize(x):
    return (x-np.min(x))/(np.max(x)-np.min(x))

# Make enthalpy temperature plot


def plot_HT(txtfile, label, color):
    data = np.loadtxt(txtfile, skiprows=1)
    temperature = data[:, 0]
    enthalpy = normalize(data[:, 1])
    plt.plot(temperature, enthalpy, color, label=label)
    plt.xlabel('temperature [K]')
    plt.ylabel('Normalized enthalpy [constant * keV]')
    plt.title('Enthalpy of Al at different temperatures, with the Al_mm potential and presimulation: NVT thermalisation')


def plot_HT2(T, H, label, color):

    temperature = T
    enthalpy = normalize(H)
    plt.plot(temperature, enthalpy, color, label=label)
    plt.xlabel('temperature [K]')
    plt.ylabel('Normalized enthalpy [constant * keV]')
    plt.title('Enthalpy Al, Al_mm potential, with NVT thermalisation')


# plot_HT('HT_data_444.txt', 'Al 4x4x4', 'ro')
# plot_HT('HT_data_777.txt', 'Al 7x7x7', 'go')
# plot_HT('HT_data_101010.txt', 'Al 10x10x10', 'bo')

# plt.legend()
# plt.show()

# Average based on decimal cut-off
def averaging_system(x, y):
    x_uniques, indices = np.unique(
        np.around(x, decimals=0), return_inverse=True)
    AV = np.zeros(len(x_uniques))
    for index, element in enumerate(x_uniques):
        temporary = []
        for indexHT, T in enumerate(x):
            if (T <= element + 0.5) and (T >= element - 0.5):
                temporary.append(y[indexHT])
        AV[index] = np.average(temporary)
    return x_uniques, AV

# Average out enthalpy every 1K


def averaged_enthalpy(txtfile, averagingfactor=1):
    data = np.loadtxt(txtfile, skiprows=1)
    temperature = data[:, 0]/averagingfactor
    enthalpy = data[:, 1]
    x, y = averaging_system(temperature, enthalpy)
    return x*averagingfactor, y


# plot the derivative of the enthalpy
def plot_dHdT(T, H, label, color):
    temperature = T
    enthalpy = H
    dHdT = np.gradient(enthalpy, temperature)
    plt.plot(temperature, dHdT, color, label=label)
    plt.xlabel('temperature [K]')
    plt.ylabel('dH/dT [keV/K]')
    plt.title('Derivative of the enthalpy of Al at different temperatures, with the Al_mm potential and presimulation: NVT thermalisation')


def determine_melting_point_and_error(T: np.array, H: np.array, stdv_param: float = 3, start_outlier: int = 1, end_outlier: int = 1):
    # take the derivative of the enthalpy
    dHdT = np.gradient(H, T)

    # determine the average and standard deviation of the derivative for selected temperature ranges (1/3 left and 1/3 right)
    Lavg, Lstdv = np.average(
        dHdT[:round(len(dHdT)/3)]), np.std(dHdT[:round(len(dHdT)/3)])
    Ravg, Rstdv = np.average(
        dHdT[2*round(len(dHdT)/3):]), np.std(dHdT[2*round(len(dHdT)/3):])

    # check which points lay beyond the standarddeviation
    L_outliers_T = np.empty(0)
    R_outliers_T = np.empty(0)
    L_outliers_H = np.empty(0)
    R_outliers_H = np.empty(0)
    for index, point in enumerate(list(zip(T, dHdT))):
        # determine left outliers
        if (point[1] > Lavg + stdv_param*Lstdv) or (point[1] < Lavg - stdv_param*Lstdv):
            # Only take into account the middle part of the data
            if index > round(len(dHdT)/3) and index < 2*round(len(dHdT)/3):
                L_outliers_T = np.append(L_outliers_T, point[0])
                L_outliers_H = np.append(L_outliers_H, point[1])
        # determine right outliers
        if (point[1] > Ravg + stdv_param*Rstdv) or (point[1] < Ravg - stdv_param*Rstdv):
            # Only take into account the middle part of the data
            if index > round(len(dHdT)/3) and index < 2*round(len(dHdT)/3):
                R_outliers_T = np.append(R_outliers_T, point[0])
                R_outliers_H = np.append(R_outliers_H, point[1])

    # determine melting point and error
    assert len(
        L_outliers_T) >= start_outlier, f"start_outlier = {start_outlier} is too large, there are only {len(L_outliers_T)} outliers, lower selection criterion by lowering stdv_param"
    first_outlier = heapq.nsmallest(start_outlier, L_outliers_T)[-1]
    assert len(
        R_outliers_T) >= end_outlier, f"end_outlier = {end_outlier} is too large, there are only {len(R_outliers_T)} outliers, lower selection criterion by lowering stdv_param"
    last_outlier = heapq.nlargest(end_outlier, R_outliers_T)[-1]
    melting_point = (first_outlier + last_outlier)/2
    error = (last_outlier - first_outlier)/2
    # plt.show()
    # plt.plot(L_outliers_T, L_outliers_H, 'yo')
    # plt.plot(R_outliers_T, R_outliers_H, 'ro')
    # plt.show()
    return melting_point, error

# Give averaged data, MP and error


def results(T, H, stdv_param=3, start_outlier=1, end_outlier=1):
    plot_HT2(T, H, 'Processed data 4x4x4', 'bo')
    MP, ER = determine_melting_point_and_error(
        T, H, stdv_param, start_outlier, end_outlier)
    plt.errorbar(MP, 0.5, xerr=ER, fmt='o', color='black', label='Melting point: ' +
                 str(round(MP, 2)) + ' +/- ' + str(round(ER, 2)) + ' K')
    return MP, ER
    # print the melting point and error in the legend
    # plt.legend(['Al 10x10x10', 'Melting point: ' + str(round(MP, 2)) + ' +/- ' + str(round(ER, 2)) + ' K'])

# # plot_dHdT(*averaged_enthalpy('HT_data_101010.txt', 20), 'Al 10x10x10', 'mo')
# plot_HT2(*averaged_enthalpy('HT_data_101010.txt',20), 'Al 10x10x10', 'bo')
# MP, ER = determine_melting_point_and_error(
#     *averaged_enthalpy('HT_data_101010.txt', 20), stdv_param=3, start_outlier=2, end_outlier=2)


# results(*averaged_enthalpy(txtfile='HT_data_101010.txt', averagingfactor=20), 3, 2, 2)
# plt.show()

########################################### HEATING WITH SMOOTHING: Savitzky-Golay filter   ###################################################################

def Savitzky_Golay_filter(T, H):
    # Apply moving average to smooth the data
    window_size = 3
    smooth_enthalpy = np.convolve(H, np.ones(
        window_size)/window_size, mode='valid')
    smooth_temperature = T[window_size//2:-(window_size//2)]
    return smooth_temperature, smooth_enthalpy


def Savitzky_Golay_filter2(T, H):
    H_smooth = signal.savgol_filter(
        H, window_length=11, polyorder=3, mode="nearest")
    return T, H_smooth


def gradient(T, H):
    return T, np.gradient(H, T)

# plot_HT2(*averaged_enthalpy('HT_data_777.txt', 10), 'Al 7x7x7', 'ro')
# plot_HT2(*gradient(*Savitzky_Golay_filter(*averaged_enthalpy('HT_data_777.txt', 10))), 'Al 7x7x7', 'yo')
# plot_HT2(*gradient(*Savitzky_Golay_filter(*Savitzky_Golay_filter(*averaged_enthalpy('HT_data_777.txt', 10)))), 'Al 7x7x7', 'bo')
# plot_HT2(*gradient(*Savitzky_Golay_filter(*Savitzky_Golay_filter(*Savitzky_Golay_filter(*averaged_enthalpy('HT_data_777.txt', 10))))), 'Al 7x7x7', 'ko')

# plot_HT2(*gradient(*Savitzky_Golay_filter2(*Savitzky_Golay_filter(*averaged_enthalpy('HT_data_777.txt', 10)))), 'Al 7x7x7', 'yo')
# plot_HT2(*gradient(*Savitzky_Golay_filter2(*Savitzky_Golay_filter2(*Savitzky_Golay_filter2(*averaged_enthalpy('HT_data_777.txt', 10))))), 'Al 7x7x7', 'ro')

# 2nd filter is the best one

# plt.show()


# PPT plot:


def filter_vs_no_filter():
    T, H = averaged_enthalpy('HT_data_101010.txt', 10)
    plt.plot(*gradient(T, H), 'yo', label='Unfiltered data')
    plt.plot(*gradient(*Savitzky_Golay_filter2(*
             Savitzky_Golay_filter2(T, H))), 'bo', label='Filtered data')
    plt.legend()
    plt.xlabel('temperature [K]')
    plt.ylabel('dH/dT [keV/K]')
    plt.show()
    plt.title('Enthalpy gradient filter comparison')
# filter_vs_no_filter()


################################################ FINAL RESULTS #####################################################################

def Normalize(T, H):
    return T, normalize(H)


def final_results1(txtfile, title):
    # Average data
    Tav, Hav = averaged_enthalpy(txtfile, 10)
    # Filter data
    T, H = Savitzky_Golay_filter2(
        *Savitzky_Golay_filter2(Tav, Hav))
    # Uncomment for gradient
    # plt.plot(*Normalize(*gradient(T, H)), 'ro', label='gradient')

    # Plot processed data
    plt.plot(T, normalize(H), 'bo', label='Processed data' +
             ' (' + str(title) + ')')
    MP, ER = determine_melting_point_and_error(T, H)
    # Plot melting point and error in the middle of the plot
    plt.errorbar(MP, 0.5, xerr=ER, fmt='o', color='black', label='Melting point: ' +
                 str(round(MP, 2)) + ' +/- ' + str(round(ER, 2)) + ' K' ' (' + str(title) + ')')
    # Plot experimental melting point
    # plt.axvline(x=933, color='black', linestyle='--', label='Experimental melting point of Al')
    # plt.title(title)
    # plt.xlabel('temperature [K]')
    # plt.ylabel('Normalized enthalpy')
    # plt.legend()


def final_results2(txtfile, title):
    # Average data
    Tav, Hav = averaged_enthalpy(txtfile, 10)
    # Filter data
    T, H = Savitzky_Golay_filter2(
        *Savitzky_Golay_filter2(Tav, Hav))
    # Uncomment for gradient
    # plt.plot(*Normalize(*gradient(T, H)), 'ro', label='gradient')

    # Plot processed data
    plt.plot(T, normalize(H), 'ro', label='Processed data' +
             ' (' + str(title) + ')')
    MP, ER = determine_melting_point_and_error(T, H)
    # Plot melting point and error in the middle of the plot
    plt.errorbar(MP, 0.6, xerr=ER, fmt='o', color='green', label='Melting point: ' +
                 str(round(MP, 2)) + ' +/- ' + str(round(ER, 2)) + ' K' ' (' + str(title) + ')')
    # plt.title(title)
    plt.xlabel('temperature [K]')
    plt.ylabel('Normalized enthalpy')
    # plt.legend()

# Plot experimental melting point

def nice_plot():
    final_results1('HT_data_101010.txt', '10x10x10')
    final_results2('HT_data_777.txt', '7x7x7')
    plt.axvline(x=933, color='black', linestyle='--',
                label='Experimental melting point: 933K')
    plt.legend()
    plt.show()

def linear_regression(x, y):
    a, b = np.polyfit(x, y, deg=1)
    return a, b
