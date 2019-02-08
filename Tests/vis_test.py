import sys
sys.path.append('../Rocket/')
sys.path.append('../Forces/')
sys.path.append('../Trajectory/')
sys.path.append('../Visual/')
import numpy as np
import Forces as Forces
import Trajectory as Trajectory
import matplotlib.pyplot as plt
import visual

from Rocket1 import RocketSimple

rad2deg = 180/np.pi
deg2rad = np.pi/180

def run(int_inclination, ramp_length, time_step, sim_time):
    rocket_file = 'myRocket.dot'
    path = 'myRocket1/'
    rocket1 = RocketSimple.from_file(rocket_file, path)
    int_inclination = 10/180.0*np.pi
    ramp_length = 2.5*rocket1.getLength()
    time_step = 0.005
    sim_time = 15


    (t, position, euler, linearVelocity, angularVelocity, AoA, thrust, gravity, drag, lift) \
    = Trajectory.calculateTrajectory(rocket1, int_inclination, ramp_length, time_step, sim_time)

    sample_rate = 1 / time_step

    t = 0
    COM = []
    COP = []
    time_domain = []
    while t < sim_time:
        time_domain = time_domain + [t]
        t += time_step

    for t in time_domain:
        COM = COM + [rocket1.getCOM(t)[0]]
    for i in range(len(AoA)):
        #COP = COP + [rocket1.getCOP(AoA[i], np.linalg.norm(linearVelocity[i]))[0]]
        COP = COP + [0]


    visual.launch(sample_rate, position, euler, COM, COP, thrust, gravity, lift, drag)

run(6, 520, 0.005, 60)
