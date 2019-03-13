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
simTime = 45 # [sec]

windObj = Wind.pinkWind(29.5, 4, np.pi, simTime + 0.1)

# Calculating trajectory
tra1 = Trajectory.calculateTrajectory(
rocketObj, initialInclination, rampLength, timeStep, simTime)

tra2 = Trajectory.calculateTrajectory(
rocketObj, initialInclination, rampLength, timeStep, simTime, wind = windObj)

t = tra1[0]

x1 = tra1[1][:,0]
y1 = tra1[1][:,1]
z1 = -tra1[1][:,2]

aoa1 = tra1[3]

w_speed0 = np.sqrt(tra1[10][:,0]**2 + tra1[10][:,1]**2 + tra1[10][:,2]**2) # Just for debugging

x2 = tra2[1][:,0]
y2 = tra2[1][:,1]
z2 = -tra2[1][:,2]

aoa2 = tra2[3]

w_speed = np.sqrt(tra2[10][:,0]**2 + tra2[10][:,1]**2 + tra2[10][:,2]**2)

plt.plot(np.sqrt(x1**2 + y1**2), z1, label="Trajectory with no wind")
plt.plot(np.sqrt(x2**2 + y2**2), z2, label="Trajectory affected by wind")
plt.title("Trajectory")
plt.xlabel("Lateral distance from launch ramp [m]")
plt.ylabel("Altitude [m]")
plt.grid()
plt.legend()
#plt.gcf().gca().set_aspect(1)
plt.show()

plt.plot(t, w_speed0)
plt.plot(t, w_speed)
plt.show()

plt.plot(t, aoa1)
plt.plot(t, aoa2)
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
