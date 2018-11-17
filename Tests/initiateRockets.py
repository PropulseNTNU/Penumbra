import sys
sys.path.append('../Rocket/')
from Rocket2 import Rocket
from Rocket1 import RocketSimple

# FOR ROCKET CLASS 1
rocket_file = 'myRocket.dot'
path = 'myRocket1/'
myRocket1 = RocketSimple.from_file(rocket_file, path)
# myRocket1.printSpecifications(0)

# FOR ROCKET CLASS 2
sample_file = 'full_report_edited.out'
init_file = 'initFile.dot'
path = 'myRocket2/'
myRocket2 = Rocket.from_file_without_AoAspeed(init_file, sample_file, path)
myRocket2.plot()
