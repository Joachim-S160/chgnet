"""

Description:
Get data from trajectory files
and return physical properties

"""

import numpy as np
import matplotlib.pyplot as plt
from chgnet.model import CHGNet
from ase.io.trajectory import Trajectory
from pymatgen.io.ase import AseAtomsAdaptor



def normalize(x: np.ndarray):
    """
    Args: array x
    returns: normalized array x    
    """
    return (x-np.min(x))/(np.max(x)-np.min(x))

def get_temperature(file: str):
    """
    input = trajectory file
    output = 1D array of temperature in Kelvin
    """
    traj = Trajectory(file)[:]
    return np.array([atoms.get_temperature() for atoms in traj])

def get_time(file: str):
    """
    input = trajectory file
    output = 1D array of time in picoseconds
    """
    traj = Trajectory(file)[:]
    # * 0.2 to convert from fs to ps (timesep is 2 fs and trajectory saves data every loginterval = 100 steps) so * 2 * 100 /1000 to get ps
    return np.arange(len(traj)) * 0.2

def get_volume(file: str):
    """
    input = trajectory file
    output = 1D array of volume in A^3
    """
    traj = Trajectory(file)[:]
    return np.array([atoms.get_volume() for atoms in traj])

def get_density(file: str):
    """
    input = trajectory file
    output = 1D array of density in g/cm^3
    """
    traj = Trajectory(file)[:]
    Volume = get_volume(file)
    mass = traj[0].get_masses().sum()
    # Convert density from amu/A^3 to g/cm^3 by multiplying by 1.66054
    return 1.66054 * mass / Volume

def get_pressure(file: str):
    """
    input = trajectory file,
    output = 1D array of pressure in GPa
    """
    # stress tensors = sigma_ij at each frame
    # sigma_ij = sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_xz, sigma_yz
    # pressure = - 1/3 * (sigma_xx + sigma_yy + sigma_zz)
    traj = Trajectory(file)[:]
    stresstensors = [atoms.get_stress(include_ideal_gas=True) for atoms in traj]
    pressure = np.array([-1/3 * (sigma[0] + sigma[1] + sigma[2]) for sigma in stresstensors]) 
    # Conversion from eV/A^3 to GPa:
    pressure = pressure * 160.21766208
    return pressure

def get_total_energy(file: str):
    """
    input: trajectory file
    output: 1D array of total energy in eV
    """
    traj = Trajectory(file)[:]
    return np.array([atoms.get_total_energy() for atoms in traj])

def get_enthalpy(file: str, one_atm_pressure:bool=True):
    """
    input = trajectory file,
            one_atm_pressure = if True pressure is 1 atm, if False the real variable pressure is used 
    output = 1D array of enthalpy in eV
    """
    if one_atm_pressure:
        # Convert 1 atm to eV/A^3: atm to GPA is 10**(-4) and GPa to eV/A^3 is 1/160
        pressure = 1/160 * 10**(-4) * np.ones_like(get_pressure(file))
    else:
        pressure = get_pressure(file)
    Etot = get_total_energy(file)
    Volume = get_volume(file)
    return normalize(Etot + pressure * Volume)

def get_kinetic_energy(file: str):
    """
    input = trajectory file
    output = 1D array of kinetic energy in eV
    """
    return np.array([atoms.get_kinetic_energy() for atoms in Trajectory(file)[:]])

def get_potential_energy(file: str):
    """
    input = trajectory file
    output = 1D array of potential energy in eV
    """
    return np.array([atoms.get_potential_energy() for atoms in Trajectory(file)[:]])

def time_averaging(quantity: np.ndarray, window_factor: float = 0.5, stdv: bool = False):
    """
    Args:   
            1D array of quantity to be averaged
            window_factor = factor to determine the size of the window
    returns:
            averaged quantity 
    """
    assert window_factor <= 1, "window_factor must be less than 1"
    if stdv:
        return np.std(quantity[int(np.round(len(quantity)*window_factor)):])
    return np.average(quantity[int(np.round(len(quantity)*window_factor)):])

def convergence_error(quantity: np.ndarray, window_factor: float = 0.5):
    """
    Args:   
            1D array of quantity to be averaged
            window_factor = factor to determine the size of the window
    returns:
            convergence error of the quantity: (min and max)/2 in the window
    """
    assert window_factor <= 1, "window_factor must be less than 1"
    window = quantity[int(np.round(len(quantity)*window_factor)):]
    return (np.max(window) - np.min(window))/2