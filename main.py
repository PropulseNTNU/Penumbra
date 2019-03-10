"""
Script for running one complete simulation
    - Initialize a rocket
    - Calculate trajectory
    - Visualize trajectory
    - Plot results

Version: WIP
Last edit: 08.02.2019

--Propulse NTNU--
"""

# TODO
# Add visual part

import sys
sys.path.append('Rocket/')
sys.path.append('Forces/')
sys.path.append('Trajectory/')
sys.path.append('Visual')
import numpy as np
from Rocket1 import RocketSimple
from Rocket2 import Rocket
import Trajectory
import Kinematics
import matplotlib.pyplot as plt
#import visual

rad2deg = 180/np.pi

# Initialize a rocket
# FOR ROCKET CLASS 1
path1 = 'Tests/myRocket1/'
init_file1 = 'myRocket.dot'
Rocket1 = RocketSimple.from_file(init_file1, path1)

# Specify initial conditions
initialInclination = 10/180.0*np.pi
launchRampLength = 2.0*Rocket1.getLength()
timeStep = 0.05
simulationTime= 25
trajectory = Trajectory.calculateTrajectoryWithBrakes(Rocket1, initialInclination, launchRampLength,
                                            timeStep, simulationTime)
# Kinematics
t = trajectory[0]
position = trajectory[1]
orientation = trajectory[2]
AoA = trajectory[3]
linearVelocity = trajectory[4]
angularVelocity = trajectory[5]

#Forces (as len(t)x3 matrices, one row correspond to 1 instance in time )
drag = trajectory[6]
lift = trajectory[7]
gravity = trajectory[8]
thrust = trajectory[9]

# Visualize trajectory

# Plot
# Position
plt.figure()
ax1 = plt.subplot(311, ylabel='x [m]')
ax1.plot(t, position[:,0], label='x', lw=2, c='r')
ax1.grid()
plt.subplots_adjust(hspace=0.5)
ax2 = plt.subplot(312, ylabel='y [m]')
ax2.plot(t, position[:,1], label='y', lw=2, c='b')
ax2.grid()
plt.subplots_adjust(hspace=0.5)
ax3 = plt.subplot(313, xlabel='time [s]', ylabel=' altitude [m]')
ax3.plot(t, -position[:,2], label='h', lw=2, c='g')
ax3.grid()
plt.figure()
ax1 = plt.subplot(111, xlabel='r [m]', ylabel=' h [m]')
ax1.plot(np.sqrt(position[:,0]**2+position[:,1]**2), -position[:,2], lw=2, c='r')
ax1.set_title('h vs. r')
ax1.grid()
ax1.axis('equal')

# Orientation
plt.figure()
plt.title('Orientation')
ax1 = plt.subplot(411, ylabel='pitch [deg]')
ax1.plot(t, orientation[:,0]*rad2deg, label='pitch', lw=2, c='r')
ax1.grid()
plt.subplots_adjust(hspace=0.5)
ax2 = plt.subplot(412, ylabel='yaw [deg]')
ax2.plot(t, orientation[:,1]*rad2deg, label='yaw', lw=2, c='b')
ax2.grid()
plt.subplots_adjust(hspace=0.5)
ax3 = plt.subplot(413, ylabel='roll [deg]')
ax3.plot(t, orientation[:,2]*rad2deg, label='roll', lw=2, c='g')
ax3.grid()
plt.subplots_adjust(hspace=0.5)
ax1 = plt.subplot(414, xlabel='time [s]', ylabel='AoA [deg]')
ax1.plot(t[15:], AoA[15:]*rad2deg, label='AoA', lw=2, c='k')
ax1.grid()

# Linear Velocity
plt.figure()
ax1 = plt.subplot(311, ylabel='v_x [m/s]')
ax1.plot(t, linearVelocity[:,0], lw=2, c='r')
ax1.grid()
ax1.set_title('Velocity (world coords)')
plt.subplots_adjust(hspace=0.5)
ax2 = plt.subplot(312, ylabel='v_y [m/s]')
ax2.plot(t, linearVelocity[:,1], lw=2, c='b')
ax2.grid()
plt.subplots_adjust(hspace=0.5)
ax3 = plt.subplot(313, xlabel='time [s]', ylabel='v_altitude [m/s]')
ax3.plot(t, -linearVelocity[:,2], lw=2, c='g')
ax3.grid()

# Angular velocity
plt.figure()
ax1 = plt.subplot(311, ylabel='w_x [rad/s]')
ax1.plot(t, angularVelocity[:,0], label='roll rate', lw=2, c='r')
ax1.grid()
ax1.legend(loc='upper right')
plt.subplots_adjust(hspace=0.5)
ax2 = plt.subplot(312, ylabel='w_y [rad/s]')
ax2.plot(t, angularVelocity[:,1], label='pitch rate', lw=2, c='b')
ax2.grid()
ax2.legend(loc='upper right')
plt.subplots_adjust(hspace=0.5)
ax3 = plt.subplot(313, xlabel='time [s]', ylabel='w_z [rad/s]')
ax3.plot(t, angularVelocity[:,2], label='yaw rate', lw=2, c='g')
ax3.grid()
ax3.legend(loc='upper right')

# Forces
plt.figure()
ax1 = plt.subplot(311, ylabel='[N]')
ax1.plot(t, thrust[:,0], label='thrust-x', lw=2, c='r')
ax1.plot(t, drag[:,0], label='drag-x', lw=2, c='b')
ax1.plot(t, lift[:,0], label='lift-x', lw=2, c='g')
ax1.plot(t, gravity[:,0], label='gravity-x', lw=2, c='k')
ax1.set_title('Forces (body-coords)')
ax1.legend(loc='upper right')
ax1.grid()
plt.subplots_adjust(hspace=0.5)
ax2 = plt.subplot(312, ylabel='[N]')
ax2.plot(t, thrust[:,1], label='thrust-y', lw=2, c='r')
ax2.plot(t, drag[:,1], label='drag-y', lw=2, c='b')
ax2.plot(t, lift[:,1], label='lift-y', lw=2, c='g')
ax2.plot(t, gravity[:,1], label='gravity-y', lw=2, c='k')
ax2.grid()
ax2.legend(loc='upper right')
plt.subplots_adjust(hspace=0.5)
ax3 = plt.subplot(313, xlabel='time [s]', ylabel='[N]')
ax3.plot(t, thrust[:,2], label='thrust-z', lw=2, c='r')
ax3.plot(t, drag[:,2], label='drag-z', lw=2, c='b')
ax3.plot(t, lift[:,2], label='lift-z', lw=2, c='g')
ax3.plot(t, gravity[:,2], label='gravity-z', lw=2, c='k')
ax3.legend(loc='upper right')
ax3.grid()

# Show plots
plt.show()
