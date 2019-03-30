"""
Script for running one complete simulation
    - Initialize a rocket
    - Calculate trajectory
    - Visualize trajectory
    - Plot results

Version: WIP
Last edit: 12.03.2019

--Propulse NTNU--
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
sys.path.append('Rocket/')
sys.path.append('Forces/')
sys.path.append('Trajectory/')
sys.path.append('Visual')
from Rocket1 import RocketSimple
import Trajectory
import TrajectoryWithBrakes

rad2deg = 180/np.pi
deg2rad = 1/rad2deg

# Initialize a rocket
# FOR ROCKET CLASS 1
path = 'rockets/'
init_file = 'init_rocket.r'
rocket = RocketSimple.from_file(init_file, path)

# Specify initial conditions

initialInclination = 10/180.0*np.pi
launchRampLength = 5
timeStep = 0.03
simulationTime= 30

tol = 0.2  # how close to targetApogee is acceptable
targetApogee = 3048 # desired apogee
rho = 1.225 
C_d = 1.28
m = rocket.getMass(6.9)
A_max = 0.006636
C_Max = (1/(2))*rho*C_d*A_max
Cbrakes = 0.4 * C_Max # force coefficient on brakes (i.e. F_brakes = [-Cbrakes*airSpeed**2,0,0])
print("CBrakes: ", Cbrakes)

trajectory = TrajectoryWithBrakes.calculateTrajectoryWithBrakes(rocket, initialInclination, launchRampLength,timeStep, simulationTime, tol, targetApogee, Cbrakes)

#initialInclination = 0*deg2rad
#launchRampLength = 5
#timeStep = 0.03
#simulationTime = 30
#trajectory = Trajectory.calculateTrajectory(rocket, initialInclination,launchRampLength,timeStep, simulationTime)

# Kinematics
t = trajectory[0]
position = trajectory[1]
orientation = trajectory[2]
AoA = trajectory[3]
linearVelocity = trajectory[4] #for some reason this is velocity in world frame?
angularVelocity = trajectory[5]

lookUpTable=[]
hoydeN=0
hoyde=0
diff=0

for i in range(len(position[:,2])):
   hoydeN=-int(np.floor(position[:,2][i]))
   print("Hoyde: ", hoyde)
   print("HoydeN: ", hoydeN)
   if hoydeN < hoyde:
       continue
   if hoydeN>hoyde:
       #var2=int(np.floor(position[:,2][i+1]))
       diff=hoydeN-hoyde
       print("Diff: ", diff)
       if i==0:
           a=(-linearVelocity[:,2][i])/diff
           for j in range(diff):
               lookUpTable.append(a*j)
               print("added: ", a*j)
               hoyde+=1
       else:
           a=(-linearVelocity[:,2][i]-(-linearVelocity[:,2][i-1]))/diff
           for j in range(diff):
               lookUpTable.append(a*j+lookUpTable[hoyde-j-1])
               print("added: ", a*j + lookUpTable[hoyde -j -1])
               hoyde+=1
   lookUpTable.append(-linearVelocity[:,2][i])
   print("Added: ", -linearVelocity[:,2][i])
   hoyde+=1
print(len(lookUpTable))
print(lookUpTable)

#Forces (as len(t)x3 matrices, one row correspond to 1 instance in time )
drag = trajectory[6]
lift = trajectory[7]
gravity = trajectory[8]
thrust = trajectory[9]
#aero_coeff = trajectory[11]

# Print trajectory statistics
Trajectory.printTrajectoryStatistics(rocket, position, linearVelocity, t)
# Plot
# Force coefficients Cd and Cn
plt.plot(-linearVelocity[:,2]/343, -drag[:,0], label='Drag(t)', lw=2, c='r')
#plt.plot(t, aero_coeff[:,1], label='Cn(t)', lw=2, c='b')
plt.grid()
#plt.ylim(0.2, 1)
plt.legend()
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
ax3 = plt.subplot(313, xlabel='time [s]', ylabel='v_z [m/s]')
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
 
#Show plots
plt.show()