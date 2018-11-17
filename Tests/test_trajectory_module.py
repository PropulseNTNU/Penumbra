import sys
sys.path.append('../Rocket/')
sys.path.append('../Forces/')
sys.path.append('../Trajectory/')
import Forces as Forces
import Trajectory as Trajectory
import Rocket2 as Rocket

def test_trajectory_module():
    sample_file = 'full-report.out'
    init_file = 'initFile.dot'
    path = 'myRocket2/'
    rocket = Rocket.from_file_with_AoAspeed(init_file, sample_file, path)
    initialInclination = 6.0/180.0*np.pi
    launchRampLength = 3*rocket.getHeight()
    timeStep = 0.1
    simulationTime= 10
    (t, position, quaternion, linearVelocity, angularVelocity, AoA, thrust, gravity, drag, lift) = Trajectory.calculateTrajectory(rocket, initialInclination, launchRampLength, timeStep, simulationTime)
    print(t)
    print(position)

def main():
    test_trajectory_module()

main()
