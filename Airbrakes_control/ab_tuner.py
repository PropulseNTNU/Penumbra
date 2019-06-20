import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('../Rocket/')
sys.path.append('../Forces/')
sys.path.append('../Trajectory/')
sys.path.append('../Optimization/')
from Rocket2 import RocketCFD
from Rocket1 import RocketSimple
import TrajectoryWithBrakes as traj
import Wind
import Optimization as optim
from airbrakes import airbrakes_main
import Kinematics

deg2rad = np.pi/180.0
rad2deg = 1/deg2rad

# Initialize a rocket (CFD)
path2 = '../rockets/sleipner_CFD/'
path1 = '../rockets/sleipner_analytic/'
drag_report = 'drag_report.txt'
lift_report = 'lift_report.txt'
moment_report = 'moment_report.txt'
init_file = 'init_rocket.r'

rocket1 = RocketSimple.from_file(init_file, path1)
rocket1.setMass(28700e-3)
rocket1.setCOM(-1270e-3)

rocket = RocketCFD.from_file(init_file, drag_report, lift_report,\
moment_report, rocket1, path2)

T = rocket.getMotor().getBurnTime()

# Specify initial conditions
initialInclination = 1e-7
rampLength = 5.2
timeStep = 0.03
simTime = 30
windObj = Wind.pinkWind(simTime, [1.5, 1.9], alt0 = 1.5, intensity = 0.1)

params = [initialInclination, rampLength, timeStep, simTime]

alt = np.zeros(int(simTime / timeStep) + 2)
acc = np.zeros(int(simTime / timeStep) + 2)

def control_function(t, r, quat, rdot, rdotdot, Cbrakes_in, Tbrakes, lookuptable):
    global T
    Rbody2Inertial = Kinematics.Rquaternion(quat)
    rdotdot_world = Rbody2Inertial @ rdotdot[0:3].T
    if int((t / 0.03)) < len(alt): alt[int(t / 0.03)] = -r[2]
    if int((t / 0.03)) < len(acc): acc[int(t / 0.03)] = np.linalg.norm(rdotdot_world)
    a = (airbrakes_main(-r[2], -rdotdot_world[2], timeStep))
    if t > T:
        print(round(t, 1), Cbrakes_in*(a / 100), sep='\t')
        return Cbrakes_in*(a / 100)
    else:
        print(round(t, 1), 0, sep='\t')
        return 0

def main():
    C = 0.001554391452
    windObj = Wind.pinkWind(simTime, [1.5,  1.9], alt0 = 1.5, intensity = 8)
    t, position, euler, AoA, velocity, angularVelocity, drag, lift, gravity,thrust, windVelocities =\
    traj.calculateTrajectoryWithAirbrakes(rocket, initialInclination,\
    rampLength, timeStep, simTime, Tbrakes = T, Cbrakes_in = C, conFunc = control_function)
    fig, ax = plt.subplots(nrows=2, ncols=2)
    ax[0][0].plot(t, -position.T[2])
    ax[1][0].plot(t, alt)
    ax[0][1].plot(acc.T)
    plt.show()

main()
