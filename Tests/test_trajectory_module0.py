import sys
sys.path.append('../Rocket/')
sys.path.append('../Forces/')
sys.path.append('../Trajectory/')
import numpy as np
import Forces as Forces
import Trajectory as Trajectory
from Rocket2 import Rocket
from Rocket1 import RocketSimple
import matplotlib.pyplot as plt

def test_trajectory_module():
    rad2deg = 180/np.pi
    # Simple rocket (class 1) #
    rocket_file = 'myRocket.dot'
    path = 'myRocket1/'
    rocket1 = RocketSimple.from_file(rocket_file, path)
    #                         #
    sample_file = 'V13_CFD.txt'
    init_file = 'V13_data.dot'
    path = 'V13/'
    rocket2 = Rocket.from_file_with_AoAspeed(init_file, sample_file, path)
    initialInclination = 1/180.0*np.pi
    launchRampLength = 2*rocket1.getLength()
    timeStep = 0.002
    simulationTime= 30
    #AoA, thrust, gravity, drag, lift
    (t, position, euler, linearVelocity, angularVelocity) = \
    Trajectory.calculateTrajectory(rocket1, initialInclination, launchRampLength, timeStep, simulationTime)

    for i in range(len())

    plt.plot(np.sqrt())
    plt.show()

def main():
    test_trajectory_module()

main()
