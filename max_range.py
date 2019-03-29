'''
Alternative main file
Last edit: 02.03.2019

--Propulse NTNU--
'''
# Telling the system where to look
import sys
sys.path.append("Forces/")
sys.path.append("Rocket/")
sys.path.append("Trajectory/")
sys.path.append("Visual/")
# Import of external libs
import numpy as np
import matplotlib.pyplot as plt
# Import of internal libs
import Wind
from Rocket1 import RocketSimple
import Trajectory
import Kinematics

f = {'family' : 'arial',
        'weight' : 'normal',
        'size'   : 12}

deg2rad = np.pi / 180
rad2deg = 180 / np.pi

# Initiaton of rocket object
rocketObj_path = 'Tests/myRocket1/'
rocketObj_initFile = 'myRocket.dot'
rocketObj = RocketSimple.from_file(rocketObj_initFile, rocketObj_path)

# Defining important paramters
inc45 = 45*deg2rad # [rad from straigth up]
inc85 = 5*deg2rad # [rad from straigth up]
rampLength = 5 # [M]
timeStep = 0.03 # [sec]
simTime = 90 # [sec]

# Calculating trajectory
trajectory45 = Trajectory.calculateTrajectory(
rocketObj, inc45, rampLength, timeStep, simTime, trim = True)

trajectory85 = Trajectory.calculateTrajectory(
rocketObj, inc85, rampLength, timeStep, simTime, trim = True)

x45 = trajectory45[1][:,0]
y45 = trajectory45[1][:,1]
z45 = -trajectory45[1][:,2]

x85 = trajectory85[1][:,0]
y85 = trajectory85[1][:,1]
z85 = -trajectory85[1][:,2]

max45 = np.max(np.sqrt(x45**2 + y45**2))
max85 = np.max(np.sqrt(x85**2 + y85**2))
s = "Max range\ninc = 45: " + str(round(max45, 2)) + "m" + "\ninc = 85: " + str(round(max85, 2)) + "m"

plt.plot(np.sqrt(x45**2 + y45**2), z45, label="inclination = 45")
plt.plot(np.sqrt(x85**2 + y85**2), z85, label="inclination = 85")
plt.title("Trajectory")
plt.xlabel("Lateral distance from launch ramp [m]")
plt.ylabel("Altitude [m]")
plt.grid()
plt.legend()
plt.gcf().gca().set_aspect(1)
plt.gcf().text(0.5, 0.15, s, fontdict=f)
plt.show()

# Output
#t = trajectory[0]
#position = trajectory[1]
#orientation = trajectory[2]
#AoA = trajectory[3]
#linearVelocity = trajectory[4]
#angularVelocity = trajectory[5]
#drag = trajectory[6]
#lift = trajectory[7]
#gravity = trajectory[8]
#thrust = trajectory[9]
