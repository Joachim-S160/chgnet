import matplotlib.pyplot as plt
import numpy as np
from ase.io.trajectory import Trajectory

traj = Trajectory("chgnet/MyCHGNetCode/mdNPT2_out_TiBr4.traj")[:]

stress = np.array([atoms.get_stress() for atoms in traj])
pressure = 1/3*(stress[:,0] + stress[:,1] + stress[:,2])
Volume = np.array([atoms.get_volume() for atoms in traj])
numberofatoms = len(traj[0])

density = np.zeros(len(traj))
temperature = np.zeros(len(traj))
indices = np.zeros(len(traj))
epot = np.zeros(len(traj))
for index, atoms in enumerate(traj):

    # Calculate volume
    volume = atoms.get_volume()
    potential_energy = atoms.get_potential_energy()

    # Calculate density
    epot[index] = potential_energy
    density[index] = 1 / volume
    temperature[index] = atoms.get_temperature()
    indices[index] = index * 0.2 # in ps


def time_temperature_plot():
    plt.plot(indices, temperature, 'ro')
    plt.xlabel('time (ps)')
    plt.ylabel('temperature (K)')
    a, _ = linear_fit(indices, temperature)
    plt.legend(['T = %.2f K/ps' % a])
    plt.show()

# find a in ax + b = 0
def linear_fit(x, y):
    import numpy as np
    A = np.vstack([x, np.ones(len(x))]).T
    a, b = np.linalg.lstsq(A, y, rcond=None)[0]
    return a, b
def time_density_plot():
    plt.plot(indices, density, 'ro')
    plt.xlabel('time (ps)')
    # unit of angstrom^3
    plt.ylabel('density (A^3)')
    # a, _ = linear_fit(indices, density)
    # plt.legend(['rho = %.2f g/cm^3/ps' % a])
    plt.show()
    
def temperature_density_plot():
    plt.plot(temperature, density, 'ro')
    plt.xlabel('temperature (K)')
    plt.ylabel('density (A^3)')
    a, _ = linear_fit(temperature, density)
    plt.legend(['rho = %.2f A^3/K' % a])
    plt.show()

def temperature_potential_plot():
    plt.plot(temperature, epot, 'ro')
    plt.xlabel('temperature (K)')
    plt.ylabel('potential energy (eV)')
    # a, _ = linear_fit(temperature, epot)
    # plt.legend(['Epot = %.2f eV/K' % a])
    plt.show()

def enthalpy_temperature_plot():
    pressure = 1
    Enthalpy = epot + pressure* Volume
    plt.plot(temperature, Enthalpy, 'ro')
    plt.xlabel('temperature (K)')
    plt.ylabel('enthalpy (eV)')
    plt.show()

# temperature_potential_plot()
# time_density_plot()
# enthalpy_temperature_plot()
# temperature_density_plot()
# time_temperature_plot()

def linear_regression(x, y):
    a, b = np.polyfit(x, y, deg=1)
    return a, b

print(linear_regression(indices, temperature)[0])