"""
Script for running one complete simulation
    - Initialize a rocket
    - Calculate trajectory
    - Use state vector to re-calculate forces
    - Visualize trajectory
    - Plot results

Version: WIP
Last edit: 28.01.2019

--Propulse NTNU--
"""

# TODO
# Initialize a rocket from folder
# Make an user interface

import sys
sys.path.append('Rocket/')
sys.path.append('Forces/')
sys.path.append('Trajectory/')
sys.path.append('Visual')
import numpy as np
from Rocket1 import RocketSimple
from Rocket2 import Rocket
import Trajectory
#import visual

path1 = ''
path2 = ''
init_file1 = ''
init_file2 = ''
force_report = ''

# Initialize a rocket
# FOR ROCKET CLASS 1
Rocket1 = RocketSimple.from_file(init_file1, path1)

# Specify initial conditions
initialInclination = 10/180.0*np.pi
launchRampLength = 2.5*Rocket1.getLength()
timeStep = 0.002
simulationTime= 15
(t, position, euler, linearVelocity, angularVelocity, AoA, thrust, gravity, drag, lift) \
= Trajectory.calculateTrajectory(Rocket1, initialInclination, launchRampLength, timeStep, simulationTime)
