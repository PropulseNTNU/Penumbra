import sys
sys.path.append('../Rocket/')
sys.path.append('../Forces/')
sys.path.append('../Trajectory/')
import numpy as np
import Forces as Forces
import Trajectory as Trajectory
from Rocket2 import Rocket
import matplotlib.pyplot as plt

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

    xs, ys, zs = position.T[0], position.T[1], position.T[2]
    hs = -zs
    rs = [np.sqrt(xs[i] ** 2 + ys[i] ** 2) for i in range(len(xs))]
    return hs, rs

def imp(filename):
    with open(filename, 'r') as f:
        fstream = f.read().strip().split('\n')
    fstream = [element.split(',') for element in fstream]
    hs = [float(fstream[i][1]) for i in range(len(fstream))]
    rs = [float(fstream[i][2]) for i in range(len(fstream))]
    return hs, rs


s_hs, s_rs = run(6, 520, 0.05, 60)
o_hs, o_rs = imp('test.csv')
plot2, = plt.plot(s_rs, s_hs)
plot1, = plt.plot(o_rs, o_hs)
plt.title("Trajectory (separate time domains)")
plt.grid()
plt.ylabel("Height [m]")
plt.xlabel("Lateral distance [m]")
plt.gcf().set_facecolor('lightgrey')
plt.legend([plot1, plot2], ["OpenRocket", "Native"])
plt.show()
