import numpy as np
import matplotlib.pyplot as plt
from chgnet.model import CHGNet
from ase.io.trajectory import Trajectory
from pymatgen.io.ase import AseAtomsAdaptor

"""
This script was made with the sole purpose of verifying the pressure units used in the data processing.
"""

file = "chgnet/MyCHGNetCode/mdNPT2_out_WCl6.traj"

# take 1 frame from the trajectory
frame = Trajectory(file)[0]

# First method:

def get_pressure(frame):
    # stress tensors = sigma_ij at each frame
    # sigma_ij = sigma_xx, sigma_yy, sigma_zz, sigma_xy, sigma_xz, sigma_yz
    # pressure = - 1/3 * (sigma_xx + sigma_yy + sigma_zz)
    stresstensors = frame.get_stress()
    sigma = stresstensors
    pressure = -1/3 * (sigma[0] + sigma[1] + sigma[2])
    return pressure

print("Pressure from get_pressure function [eV/A^3]: ", get_pressure(frame))
print("Pressure from get_pressure function [GPa]: ", get_pressure(frame) * 160.21766208)

def get_pressure2(frame):
    chgnet = CHGNet.load()
    atoms =frame
    stresstensors = chgnet.predict_structure(AseAtomsAdaptor.get_structure(atoms))['s'] 
    pressure = - stresstensors.trace() / 3
    return pressure
    
print("Pressure from get_pressure2 function [GPa]: ", get_pressure2(frame))
