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

chgnet = CHGNet.load()
# normalize function
def normalize(x):
    return (x-np.min(x))/(np.max(x)-np.min(x))

def get_THtrp(files: list, Fixed_pressure:bool=False):
    """
    input = list of trajectory files, 
            Fixed_pressure = True if pressure is fixed, False if pressure is variable
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
        
        
        
        stresstensors = [chgnet.predict_structure(AseAtomsAdaptor.get_structure(atoms))['s'] for atoms in traj]
        pressure = np.array([-1/3 * (sigma[0] + sigma[1] + sigma[2]) for sigma in stresstensors]) 
        # Make sure the pressure is 1 atm for the enthalpy calculation if need be
        # 1 atm = 10^-4 * 1/160 eV/A^3
        _pressure_ = (160)*10**(4)*np.ones_like(pressure) if Fixed_pressure else pressure
        Enthalpy = normalize(Etot + _pressure_ * Volume)
        mass = traj[0].get_masses().sum()
        print(index, file)

        # Fill the 3D array with data
        THtrfiles[index, 0, :num_timesteps_file] = [atoms.get_temperature() for atoms in traj]
        THtrfiles[index, 1, :num_timesteps_file] = Enthalpy
        THtrfiles[index, 2, :num_timesteps_file] = np.arange(num_timesteps_file) * 0.2
        # Convert density from amu/A^3 to g/cm^3 by multiplying by 1.66054
        THtrfiles[index, 3, :num_timesteps_file] = 1.66054 * mass / Volume
        # Convert pressure from eV/A^3 to GPa by multiplying by 160.21766208
        THtrfiles[index, 4, :num_timesteps_file] = _pressure_ * 160.21766208

    return THtrfiles

def get_tKEPE(files: list):
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