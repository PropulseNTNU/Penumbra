import sys
sys.path.append('../Rocket/')
sys.path.append('../Forces/')
sys.path.append('../Trajectory/')
sys.path.append('../Visual/')
import numpy as np
import Forces as Forces
import Trajectory as Trajectory
from Rocket2 import Rocket
import matplotlib.pyplot as plt
import visual

rad2deg = 180/np.pi
deg2rad = np.pi/180

def run(int_inclination, ramp_length, time_step, sim_time):
    sample_file = 'V13_CFD.txt'
    init_file = 'V13_data.dot'
    path = 'V13/'
    rocket = Rocket.from_file_with_AoAspeed(init_file, sample_file, path)
    initialInclination = int_inclination/180.0*np.pi

    (t, position, euler, linearVelocity, angularVelocity, AoA, thrust, gravity, drag, lift) \
    = Trajectory.calculateTrajectory(rocket, int_inclination, ramp_length, time_step, sim_time)

    sample_rate = 1 / time_step
    print(position)
    visual.launch(sample_rate, position, euler)

run(6, 520, 0.05, 60)
