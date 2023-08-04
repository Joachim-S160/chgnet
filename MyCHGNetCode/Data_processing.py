from scipy import signal
import numpy as np
import matplotlib.pyplot as plt
import heapq
from scipy.interpolate import make_interp_spline
from ase.io.trajectory import Trajectory
from plot import get_THtrp

THtrp = get_THtrp(["chgnet/MyCHGNetCode/mdNPT3_out_Al.traj"])

Temperature = THtrp[0, 0, :]

# NPT simulation so pressure is constant
Enthalpy = THtrp[0, 1, :]

plot_everything = True

################################### HEATING ######################################################################################

# normalize function
def normalize(x):
    return (x-np.min(x))/(np.max(x)-np.min(x))

# Make enthalpy temperature plot

def plot_HT(T, H, label, color):
    plt.plot(T, normalize(H), color, label=label)
    plt.xlabel('temperature [K]')
    plt.ylabel('Normalized enthalpy [constant * keV]')

if plot_everything:  
    plot_HT(Temperature, Enthalpy, 'Al', 'ro')
    plt.legend()
    plt.show()

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

# Average out enthalpy every 1K or more
def averaged_enthalpy(T, H, averagingfactor=1):
    temperature = T/averagingfactor
    enthalpy = H
    x, y = averaging_system(temperature, enthalpy)
    return x*averagingfactor, y


if plot_everything:  
    plt.plot(*averaged_enthalpy(Temperature, Enthalpy, 10), 'ro')
    plt.show()

# plot the derivative of the enthalpy
def plot_dHdT(T, H, label, color):
    temperature = T
    enthalpy = H
    dHdT = np.gradient(enthalpy, temperature)
    plt.plot(temperature, dHdT, color, label=label)
    plt.xlabel('temperature [K]')
    plt.ylabel('dH/dT [keV/K]')


if plot_everything:  
    plot_dHdT(*averaged_enthalpy(Temperature, Enthalpy, 10), 'Al', 'ro')
    plt.show()

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
    MP, ER = determine_melting_point_and_error(
        T, H, stdv_param, start_outlier, end_outlier)
    plt.errorbar(MP, 0.5, xerr=ER, fmt='o', color='black', label='Melting point: ' +
                 str(round(MP, 2)) + ' +/- ' + str(round(ER, 2)) + ' K')
    return MP, ER
    # print the melting point and error in the legend
    # plt.legend(['Al 10x10x10', 'Melting point: ' + str(round(MP, 2)) + ' +/- ' + str(round(ER, 2)) + ' K'])


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


def filter_vs_no_filter(T,H):
    T_, H_ = averaged_enthalpy(T,H, 10)
    plt.plot(*gradient(T_, H_), 'yo', label='Unfiltered data')
    plt.plot(*gradient(*Savitzky_Golay_filter2(*
             Savitzky_Golay_filter2(T_, H_))), 'bo', label='Filtered data')
    plt.legend()
    plt.xlabel('temperature [K]')
    plt.ylabel('dH/dT [keV/K]')
    plt.title('Enthalpy gradient filter comparison')
    plt.show()
    

if plot_everything:  
    filter_vs_no_filter(Temperature, Enthalpy)


################################################ FINAL RESULTS #####################################################################

def Normalize(T, H):
    return T, normalize(H)


def final_results1(T,H, title):
    # Average data
    Tav, Hav = averaged_enthalpy(T,H, 10)
    # Filter data
    T_, H_ = Savitzky_Golay_filter2(
        *Savitzky_Golay_filter2(Tav, Hav))
    # Uncomment for gradient
    # plt.plot(*Normalize(*gradient(T, H)), 'ro', label='gradient')

    # Plot processed data
    plt.plot(T_, normalize(H_), 'bo', label='Processed data' +
             ' (' + str(title) + ')')
    MP, ER = determine_melting_point_and_error(T_, H_,2)
    # Plot melting point and error in the middle of the plot
    plt.errorbar(MP, 0.5, xerr=ER, fmt='o', color='black', label='Melting point: ' +
                 str(round(MP, 2)) + ' +/- ' + str(round(ER, 1.1)) + ' K' ' (' + str(title) + ')')
    plt.title(title)
    plt.xlabel('temperature [K]')
    plt.ylabel('Normalized enthalpy')
    # plt.legend()



def nice_plot():
    final_results1(Temperature, Enthalpy,'Al')
    plt.axvline(x=933, color='black', linestyle='--',
                label='Experimental melting point: 933K')
    plt.legend()
    plt.show()

def linear_regression(x, y):
    a, b = np.polyfit(x, y, deg=1)
    return a, b


if plot_everything:  
    nice_plot()