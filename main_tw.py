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

deg2rad = np.pi / 180
rad2deg = 180 / np.pi

# Initiaton of rocket object
rocketObj_path = 'Tests/myRocket1/'
rocketObj_initFile = 'myRocket.dot'
rocketObj = RocketSimple.from_file(rocketObj_initFile, rocketObj_path)

# Defining important paramters
initialInclination = 6*deg2rad # [rad from straigth up]
rampLength = 5 # [M]
timeStep = 0.03 # [sec]
simTime = 25 # [sec]

windObj = Wind.engWind(29.5, 4, np.pi)

# Calculating trajectory
tra_noWind = Trajectory.calculateTrajectory(
rocketObj, initialInclination, rampLength, timeStep, simTime)

tra_wind = Trajectory.calculateTrajectory(
rocketObj, initialInclination, rampLength, timeStep, simTime, wind = windObj)

x = tra_noWind[1][:,0]
y = tra_noWind[1][:,1]
z = -tra_noWind[1][:,2]

xw = tra_wind[1][:,0]
yw = tra_wind[1][:,1]
zw = -tra_wind[1][:,2]

plt.plot(np.sqrt(x**2 + y**2), z)
plt.plot(np.sqrt(xw**2 + yw**2), zw)
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
