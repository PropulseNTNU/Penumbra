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
import datetime as dt
from scipy.stats.kde import gaussian_kde
import pickle
import signal
import os

verbose = True
interrupted = False


deg2rad = np.pi / 180
rad2deg = 180 / np.pi

# Initiaton of rocket object
rocketObj_path = 'Tests/myRocket1/'
rocketObj_initFile = 'myRocket.dot'
initialInclination = 6*deg2rad # [rad from straigth up]
rampLength = 5 # [M]
timeStep = 0.03 # [sec]
simTime = 40 # [sec]

rocketObj = RocketSimple.from_file(rocketObj_initFile, rocketObj_path)
params = [initialInclination, rampLength, timeStep, simTime]
windObj = Wind.pinkWind(10, 8, simTime + 0.1)

def cusp(out):
    if verbose: print(out)

class MonteCarlo:
    def __init__(self, n, rocket, params, thrustFreq = 0, thrustDev = 0, dragDev = 0, wind = Wind.nullWind()):
        self.n = n
        self.apogees = np.zeros(n)
        self.deltas = np.zeros(n)
        self.params = params
        self.rocket = rocket
        self.rocket.getMotor().setStochasticParams(thrustFreq, thrustDev)
        self.thrustFreq = thrustFreq
        self.thrustDev = thrustDev
        self.dragDev = dragDev
        self.wind = wind
        self.i = 0
        self.id = dt.datetime.now().strftime("%y%m%d%H%M%S%f")
        cusp("Monte")

    def showProgess(self):
        scaledown = 1.48
        fullbar = int(100 // scaledown)
        progress = int((self.i)/self.n * 100)
        pbar = int(progress // scaledown)
        meanApogee = np.mean(self.apogees[np.nonzero(self.apogees)])
        meanDelta = np.mean(self.deltas[np.nonzero(self.deltas)])
        timeLeft = meanDelta * (self.n - self.i)

        os.system("cls")
        print(\
        "Iterations".ljust(20) + \
        "Progress".ljust(20) + \
        "Time left".ljust(20) + \
        "Mean apogee".ljust(20))

        print('-'*71)

        print(\
        "{} / {}".format(self.i, self.n).ljust(20) + \
        "{}%".format(progress).ljust(20) + \
        "{}".format(str(dt.timedelta(seconds=int(timeLeft)))).ljust(20) + \
        "{}m".format(round(meanApogee, 2)).ljust(20))
        print('\n' + "  " + pbar*chr(9619) + (fullbar - pbar)*chr(9617) + '\n')

        print(self)


    def run(self):
        #os.system("cls")
        def signal_handler(signal, frame):
            global interrupted
            print("Process safely interrupted")
            interrupted = True

        signal.signal(signal.SIGINT, signal_handler)

        while self.i < self.n and not interrupted:
            T0 = dt.datetime.now()
            self.rocket.getMotor().refresh()
            self.wind.refresh()

            trajectory = Trajectory.calculateTrajectory(self.rocket,\
            *self.params[0:4], windObj = self.wind,\
            dragDeviation = self.dragDev)


            self.apogees[self.i] = np.max(-trajectory[1][:,2])

            with open(self.id + ".mc", "wb") as save:
                pickle.dump(self, save)

            Tf = dt.datetime.now()
            self.deltas[self.i] = (Tf - T0).total_seconds()

            self.i += 1

            self.showProgess()

        print("Process terminated, id: {}".format(self.id))

    @staticmethod
    def fromFile(id):
        with open(id + ".mc", "rb") as save:
            return pickle.load(save)

    def flush(self):
        self.i = 0
        self.apogees = np.zeros(self.n)

    def setStatus(status):
        self.status = status

    def setId(newId):
        self.id = newId

    def pangea(self):
        apogees = self.apogees[np.nonzero(self.apogees)]
        pdf = gaussian_kde(apogees)
        x = np.linspace(np.min(apogees) - 200, np.max(apogees) + 200, 1000)
        plt.plot(x, pdf(x))
        plt.fill_between(x, pdf(x), color = (0.9, 0.9, 0.9))
        plt.grid()
        plt.title("Apogee PDF")
        plt.xlabel("Apogee [m]")
        mn = np.mean(apogees)
        plt.plot(mn, 0, 'x')
        plt.show()

    def __str__(self):
        return "id:".ljust(20) + self.id + '\n' + \
        "n:".ljust(20) + str(self.n) + '\n' + \
        "Thrust deviation:".ljust(20) + str(self.thrustDev) + '\n' + \
        "Thrust frequency:".ljust(20) + str(self.thrustFreq) + '\n' + \
        "Drag deviation:".ljust(20) + str(self.dragDev) + 2*'\n' + \
        self.wind.__str__() + '\n'

#test = MonteCarlo(30000, rocketObj, params, thrustFreq = 20, thrustDev = 50,\
#dragDev = 0.016, wind = windObj)

test = MonteCarlo.fromFile("190327225234737225")
#test.flush()
test.run()
test.pangea()

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
