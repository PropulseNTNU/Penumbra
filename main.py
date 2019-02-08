"""
Script for running one complete simulation
    - Initialize a rocket
    - Calculate trajectory
    - Visualize trajectory
    - Plot results

Version: WIP
Last edit: 04.02.2019

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

path1 = 'Tests/myRocket1/'
init_file1 = 'myRocket.dot'

# Initialize a rocket
# FOR ROCKET CLASS 1
Rocket1 = RocketSimple.from_file(init_file1, path1)

# Specify initial conditions
initialInclination = 4/180.0*np.pi
launchRampLength = 2.0*Rocket1.getLength()
timeStep = 0.05
simulationTime= 30
trajectory = Trajectory.calculateTrajectory(Rocket1, initialInclination, launchRampLength, 
                                            timeStep, simulationTime)
# Kinematics
t = trajectory[0]
position = trajectory[1]
euler = trajectory[2]
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

# Orientation
plt.figure()
ax1 = plt.subplot(411, xlabel='time [s]', ylabel='pitch [deg]')
ax1.plot(t, euler[:,0]*rad2deg, label='pitch', lw=2, c='r')
ax1.set_title('pitch')
ax1.grid()
plt.subplots_adjust(hspace=0.5)
ax2 = plt.subplot(412, xlabel='time [s]', ylabel='yaw [deg]')
ax2.plot(t, euler[:,1]*rad2deg, label='yaw', lw=2, c='b')
ax2.set_title('yaw')
ax2.grid()
plt.subplots_adjust(hspace=0.5)
ax3 = plt.subplot(413, xlabel='time [s]', ylabel='roll [deg]')
ax3.plot(t, euler[:,2]*rad2deg, label='roll', lw=2, c='g')
ax3.set_title('roll')
ax3.grid()
plt.subplots_adjust(hspace=0.5)
ax1 = plt.subplot(414, xlabel='time [s]', ylabel='AoA [deg]')
ax1.plot(t[15:], AoA[15:]*rad2deg, label='AoA', lw=2, c='k')
ax1.set_title('AoA')
ax1.grid()

# Linear Velocity
plt.figure()
ax1 = plt.subplot(311, xlabel='time [s]', ylabel='v_x [m/s]')
ax1.plot(t, linearVelocity[:,0], lw=2, c='r')
ax1.set_title('v_x')
ax1.grid()
plt.subplots_adjust(hspace=0.5)
ax2 = plt.subplot(312, xlabel='time [s]', ylabel='v_y [m/s]')
ax2.plot(t, linearVelocity[:,1], lw=2, c='b')
ax2.set_title('v_y')
ax2.grid()
plt.subplots_adjust(hspace=0.5)
ax3 = plt.subplot(313, xlabel='time [s]', ylabel='v_h [m/s]')
ax3.plot(t, -linearVelocity[:,2], lw=2, c='g')
ax3.set_title('v_h')
ax3.grid()

# Angular velocity
plt.figure()
ax1 = plt.subplot(311, xlabel='time [s]', ylabel='w_x [rad/s]')
ax1.plot(t, angularVelocity[:,0], lw=2, c='r')
ax1.set_title('w_x')
ax1.grid()
plt.subplots_adjust(hspace=0.5)
ax2 = plt.subplot(312, xlabel='time [s]', ylabel='w_y [rad/s]')
ax2.plot(t, angularVelocity[:,1], lw=2, c='b')
ax2.set_title('w_y')
ax2.grid()
plt.subplots_adjust(hspace=0.5)
ax3 = plt.subplot(313, xlabel='time [s]', ylabel='w_z [rad/s]')
ax3.plot(t, angularVelocity[:,2], lw=2, c='g')
ax3.set_title('w_z')
ax3.grid()

# Forces
plt.figure()
ax1 = plt.subplot(311, xlabel='time [s]', ylabel='[N]')
ax1.plot(t, thrust[:,0], label='thrust', lw=2, c='r')
ax1.plot(t, drag[:,0], label='drag', lw=2, c='b')
ax1.plot(t, lift[:,0], label='lift', lw=2, c='g')
ax1.plot(t, gravity[:,0], label='gravity', lw=2, c='k')
ax1.set_title('Forces [x, y and z components]')
ax1.legend(loc='upper right')
ax1.grid()
plt.subplots_adjust(hspace=0.5)
ax2 = plt.subplot(312, xlabel='time [s]', ylabel='[N]')
ax2.plot(t, thrust[:,1], label='thrust', lw=2, c='r')
ax2.plot(t, drag[:,1], label='drag', lw=2, c='b')
ax2.plot(t, lift[:,1], label='lift', lw=2, c='g')
ax2.plot(t, gravity[:,1], label='gravity', lw=2, c='k')
ax2.grid()
ax2.legend(loc='upper right')
plt.subplots_adjust(hspace=0.5)
ax3 = plt.subplot(313, xlabel='time [s]', ylabel='[N]')
ax3.plot(t, thrust[:,2], label='thrust', lw=2, c='r')
ax3.plot(t, drag[:,2], label='drag', lw=2, c='b')
ax3.plot(t, lift[:,2], label='lift', lw=2, c='g')
ax3.plot(t, gravity[:,2], label='gravity', lw=2, c='k')
ax3.legend(loc='upper right')
ax3.grid()

# Show plots
plt.show()
