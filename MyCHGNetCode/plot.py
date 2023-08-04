import matplotlib.pyplot as plt
import numpy as np
from ase.io.trajectory import Trajectory


def get_THtr(files: list):
    """
    input = list of trajectory files, 
    output = 3D Matrix: 
    1st dimension is the file, 
    2nd dimension is temperature, enthalpy, time, density 
    3d dimension is the time step
    """

    # create a 3D array to store the data
    THtrfiles = np.zeros((len(files), 4))
    THtrfiles = np.zeros((len(files), 4))
    for index, file in enumerate(files):
        traj = Trajectory(file)[:]
        Volume = np.array([atoms.get_volume() for atoms in traj])
        Epot = np.array([atoms.get_potential_energy() for atoms in traj])
        # NPT simulation so pressure is constant
        Enthalpy = Epot + Volume
        print(index, file)
        THtrfiles[index, 0] = [atoms.get_temperature() for atoms in traj]
        THtrfiles[index, 1] = Enthalpy.tolist()
        THtrfiles[index, 2] = np.arrange(len(traj) * 0.2).tolist()
        THtrfiles[index, 3] = (1 / Volume).tolist()
    return THtrfiles


# traj = Trajectory("chgnet/MyCHGNetCode/mdNPT2_out_TiBr4.traj")[:]

# stress = np.array([atoms.get_stress() for atoms in traj])
# pressure = 1/3*(stress[:, 0] + stress[:, 1] + stress[:, 2])
# Volume = np.array([atoms.get_volume() for atoms in traj])
# numberofatoms = len(traj[0])

# density = np.zeros(len(traj))
# temperature = np.zeros(len(traj))
# time = np.zeros(len(traj))
# epot = np.zeros(len(traj))
# for index, atoms in enumerate(traj):

#     # Calculate volume
#     volume = atoms.get_volume()
#     potential_energy = atoms.get_potential_energy()

#     # Calculate density
#     epot[index] = potential_energy
#     density[index] = 1 / volume
#     temperature[index] = atoms.get_temperature()
#     time[index] = index * 0.2  # in ps

def time_temperature_plot(files: list):
    THtr = get_THtr(files)
    for index in range(len(files)):
        time = THtr[index, 2]
        temperature = THtr[index, 0]
        plt.plot(time, temperature, 'ro', label= index)
        plt.xlabel('time (ps)')
        plt.ylabel('temperature (K)')
        print(f"Linear regression for {files[index]} = {linear_regression(time, temperature)}")
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
  
time_temperature_plot(['chgnet/MyCHGNetCode/mdNPT2_out_HfF4.traj','chgnet/MyCHGNetCode/mdNPT2_out_LiCl.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiBr4.traj','chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj','chgnet/MyCHGNetCode/mdNPT2_out_TiI4.traj'])



# def time_density_plot():
#     plt.plot(time, density, 'ro')
#     plt.xlabel('time (ps)')
#     # unit of angstrom^3
#     plt.ylabel('density (A^3)')
#     # a, _ = linear_fit(time, density)
#     # plt.legend(['rho = %.2f g/cm^3/ps' % a])
#     plt.show()


# def temperature_density_plot():
#     plt.plot(temperature, density, 'ro')
#     plt.xlabel('temperature (K)')
#     plt.ylabel('density (A^3)')
#     a, _ = linear_fit(temperature, density)
#     plt.legend(['rho = %.2f A^3/K' % a])
#     plt.show()


# def temperature_potential_plot():
#     plt.plot(temperature, epot, 'ro')
#     plt.xlabel('temperature (K)')
#     plt.ylabel('potential energy (eV)')
#     # a, _ = linear_fit(temperature, epot)
#     # plt.legend(['Epot = %.2f eV/K' % a])
#     plt.show()


# def enthalpy_temperature_plot():
#     pressure = 1
#     Enthalpy = epot + pressure * Volume
#     plt.plot(temperature, Enthalpy, 'ro')
#     plt.xlabel('temperature (K)')
#     plt.ylabel('enthalpy (eV)')
#     plt.show()

# # temperature_potential_plot()
# # time_density_plot()
# # enthalpy_temperature_plot()
# # temperature_density_plot()
# # time_temperature_plot()




