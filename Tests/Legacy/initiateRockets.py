import sys
sys.path.append('../Rocket/')
from Rocket2 import Rocket
from Rocket1 import RocketSimple

# FOR ROCKET CLASS 1
rocket_file = 'myRocket.dot'
path = 'myRocket1/'
myRocket1 = RocketSimple.from_file(rocket_file, path)

# FOR ROCKET CLASS 2
#sample_file = 'V13_CFD.txt'
#init_file = 'V13_data.dot'
#path = 'V13/'
#rocket = Rocket.from_file_with_AoAspeed(init_file, sample_file, path)
#rocket.plot()
#rocket.getMotor().plotPerformance()
