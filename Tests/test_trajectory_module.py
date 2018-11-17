import sys
sys.path.append('../Rocket/')
sys.path.append('../Forces/')
sys.path.append('../Trajectory/')
import Forces as Forces
import Trajectory as Trajectory
from Rocket2 import Rocket
import numpy as np

def test_trajectory_module():
    sample_file = 'full_report_edited.dot'
    init_file = 'initFile.dot'
    path = 'myRocket2/'
    rocket = Rocket.from_file_with_AoAspeed(init_file, sample_file, path)
    initialInclination = 6.0/180.0*np.pi
    launchRampLength = 3*rocket.getLength()
    timeStep = 0.1
    simulationTime= 10
    (t, position, quaternion, linearVelocity, angularVelocity, AoA, thrust, gravity, drag, lift) = Trajectory.calculateTrajectory(rocket, initialInclination, launchRampLength, timeStep, simulationTime)
    print(t)
    print(position)

def main():
    test_trajectory_module()

main()
