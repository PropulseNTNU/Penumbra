"""
Script for running one complete simulation
    - Initialize a rocket
    - Calculate trajectory
    - Use state vector to re-calculate forces
    - Visualize trajectory
    - Plot results

Version: WIP
Last edit: 30.01.2019

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
import matplotlib.pyplot as plt
#import visual

path1 = 'Tests/myRocket1/'
init_file1 = 'myRocket.dot'

# Initialize a rocket
# FOR ROCKET CLASS 1
Rocket1 = RocketSimple.from_file(init_file1, path1)
# Specify initial conditions
initialInclination = 5/180.0*np.pi
launchRampLength = 2.0*Rocket1.getLength()
timeStep = 0.001
simulationTime= 2
(t, position, euler, linearVelocity, angularVelocity) \
= Trajectory.calculateTrajectory(Rocket1, initialInclination, launchRampLength, timeStep, simulationTime)

# Calculate plot data (AoA/Orientation, Forces etc.)

# Visualize trajectory

# Plot
plt.figure()
ax1 = plt.subplot(311, xlabel='time [s]', ylabel='x [m]')
ax1.plot(t, position[:,0], label='x', lw=2, c='r')
ax1.set_title('x')
ax1.grid()
plt.subplots_adjust(hspace=0.5)
ax2 = plt.subplot(312, xlabel='time [s]', ylabel='y [m]')
ax2.plot(t, position[:,1], label='y', lw=2, c='b')
ax2.set_title('y')
ax2.grid()
plt.subplots_adjust(hspace=0.5)
ax3 = plt.subplot(313, xlabel='time [s]', ylabel=' h [m]')
ax3.plot(t, -position[:,2], label='h', lw=2, c='g')
ax3.set_title('h')
ax3.grid()

plt.figure()
ax1 = plt.subplot(111, xlabel='r [m]', ylabel=' h [m]')
ax1.plot(np.sqrt(position[:,0]**2+position[:,1]**2), -position[:,2], lw=2, c='r')
ax1.set_title('h vs. r')
ax1.grid()
ax1.axis('equal')

# Show plots
plt.show()
