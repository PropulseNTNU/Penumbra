import sys
sys.path.append('../Rocket')
sys.path.append('../Forces')
sys.path.append('../Trajectory')
import Rocket.Rocket as rocket
import Forces.Forces as forces
import Trajectory.Trajectory as trajectory

def plot_rocket_quantities(rocket_file, path_to_folder=""):
    rocket = RocketSimple.from_file(rocket_file, path_to_folder)

def test_trajectory_module(rocket_file, path_to_folder=""):
    rocket = RocketSimple.from_file(rocket_file, path_to_folder)



def main():
    rocket_file = 'myRocket.dot'
    path_to_folder = 'Rocket/myRocket'
    initialInclination = 6*np.pi
    launchRampLength = 0.1
    timeStep = 0.001
    simulationTime = 2.5
    rocket = rocket.RocketSimple.from_file(rocket, path_to_folder)
    (t, position, quaternion, linearVelocity, angularVelocity, AoA, thrust, gravity, drag, lift) = calculateTrajectory(rocket, initialInclination, launchRampLength, timeStep, simulationTime)
    #plot_rocket_quantities(rocket_file)
    #test_trajectory_module(rocket_file, path_to_folder)
