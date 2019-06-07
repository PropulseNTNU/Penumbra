'''
Alternative main file
Last edit: 02.03.2019

--Propulse NTNU--
'''
# Telling the system where to look
import sys
sys.path.append("../Forces/")
sys.path.append("../Rocket/")
sys.path.append("../Trajectory/")
sys.path.append("../Visual/")
# Import of external libs
import numpy as np
import matplotlib.pyplot as plt
# Import of internal libs
import Wind
from Rocket1 import RocketSimple
import Trajectory
import TrajectoryWithBrakes as traBrakes
import Kinematics
import datetime as dt
from scipy.stats.kde import gaussian_kde
import scipy
import statistics
import pickle
import signal
import os
from scipy.interpolate import interp1d

verbose = True
interrupted = False

class ShootingOptimz:
    def __init__(self, rocket, initialInclination, rampLength, timeStep, simTime, Cbrakes_in, windObj = Wind.nullWind(), thrustFreq = 0, thrustDev = 0, dragDev = 0):
        self.rocket = rocket
        self.initialInclination = initialInclination
        self.launchRampLength = rampLength
        self.timeStep = timeStep
        self.simulationTime = simTime
        #self.Tbrakes = Tbrakes
        self.Cbrakes_in = Cbrakes_in
        self.windObj = windObj
        self.thrustFreq = thrustFreq
        self.thrustDev = thrustDev
        self.dragDev = dragDev
        self.params = [initialInclination, rampLength, timeStep, simTime]


    def run():
        pass

    def simpleShoot(self, t0, dt, targetApogee, tol):
        t = t0
        # Loop
        apogee = targetApogee + 2 * tol
        while(np.abs(targetApogee - apogee) > tol):
            apogee = np.max(-traBrakes.calculateTrajectoryWithAirbrakes(self.rocket,\
            self.initialInclination, self.launchRampLength, self.timeStep,\
            self.simulationTime, t, self.Cbrakes_in, self.dragDev,\
            self.windObj)[1][:,2])
            a = np.sign(targetApogee - apogee)
            if t != t0:
                if a != b: dt /= 2
            else:
                b = a
            t += a * dt
            print(str(t).ljust(20, ' ') + str(apogee).ljust(20, ' '))

    def stochasticShoot(self, t0, dt, targetApogee, tol, Tbrakes = 0, Cbrakes_in = 0):
        t = t0
        # Loop
        apogee = targetApogee + 2 * tol
        while(np.abs(targetApogee - apogee) > tol):
            mcResult = MonteCarlo(100, self.rocket, self.params,\
            thrustFreq = self.thrustFreq, thrustDev = self.thrustDev,\
            dragDev = self.dragDev, wind = self.windObj,\
            doSave = False, verbose = False, Tbrakes = t, Cbrakes_in = 1.0)

            mcResult.run()

            apogee = np.mean(mcResult.apogees[np.nonzero(mcResult.apogees)])
            apogee = np.mean(mcResult.apogees[0])
            #print(apogee)
            print(str(t).ljust(20, ' ') + str(apogee).ljust(20, ' '))

            a = np.sign(targetApogee - apogee)
            if t != t0:
                if a != b: dt /= 2
            else:
                b = a
            t += a * dt



class MonteCarlo:
    def __init__(self, n, rocket, params, thrustFreq = 0, thrustDev = 0,\
    dragDev = 0, wind = Wind.nullWind(), doSave = True, verbose = True,\
    Tbrakes = 60, Cbrakes_in = 0):
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
        self.doSave = doSave
        self.verbose = verbose
        self.Tbrakes = Tbrakes
        self.Cbrakes_in = Cbrakes_in
        if self.verbose: print("Monte")

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
        global interrupted
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

            trajectory = traBrakes.calculateTrajectoryWithAirbrakes(self.rocket,\
            *self.params[0:4], windObj = self.wind,\
            dragDeviation = self.dragDev, Tbrakes = self.Tbrakes,\
            Cbrakes_in = self.Cbrakes_in)


            self.apogees[self.i] = np.max(-trajectory[1][:,2])

            if self.doSave:
                self.save()

            Tf = dt.datetime.now()
            self.deltas[self.i] = (Tf - T0).total_seconds()

            self.i += 1

            if self.verbose: self.showProgess()

        interrupted = False
        if self.verbose: print("Process terminated, id: {}".format(self.id))

    @staticmethod
    def fromFile(id):
        with open(id + ".mc", "rb") as save:
            return pickle.load(save)

    def flush(self):
        self.i = 0
        self.apogees = np.zeros(self.n)

    def setStatus(self, status):
        self.status = status

    def save(self):
        with open(self.id + ".mc", "wb") as save:
            pickle.dump(self, save)

    def setId(self, newId):
        self.id = newId

    def pangea(self):
        apogees = self.apogees[np.nonzero(self.apogees)]
        mn = np.mean(apogees)
        stdDev = statistics.stdev(apogees)
        pdf = gaussian_kde(apogees)
        x = np.linspace(np.min(apogees) - 200, np.max(apogees) + 200, 1000)
        x1stdDev = np.linspace(mn - stdDev, mn + stdDev, 1000)
        x2stdDev = np.linspace(mn - 2 * stdDev, mn + 2 * stdDev, 1000)
        x3stdDev = np.linspace(mn - 3 * stdDev, mn + 3 * stdDev, 1000)
        #x4stdDev = np.linspace(mn - 4 * stdDev, mn + 4 * stdDev, 1000)
        plt.plot(x, pdf(x), color = (0, 0, 0))
        plt.fill_between(x, pdf(x), color = (0, 0, 0))
        #plt.fill_between(x4stdDev, pdf(x4stdDev), color = (0.3, 0.3, 0.3))
        stdDev3p = scipy.integrate.quad(pdf, mn - 3 * stdDev, mn + 3 * stdDev)[0]
        stdDev2p = scipy.integrate.quad(pdf, mn - 2 * stdDev, mn + 2 * stdDev)[0]
        stdDev1p = scipy.integrate.quad(pdf, mn - stdDev, mn + stdDev)[0]
        plt.fill_between(x3stdDev, pdf(x3stdDev), color = (0.5, 0.5, 0.5), label = r"$\mu \pm 3\sigma\equal {:.1f}\%$".format(stdDev3p * 100))
        plt.fill_between(x2stdDev, pdf(x2stdDev), color = (0.7, 0.7, 0.7), label = r"$\mu \pm 2\sigma\equal {:.1f}\%$".format(stdDev2p * 100))
        plt.fill_between(x1stdDev, pdf(x1stdDev), color = (0.9, 0.9, 0.9), label = r"$\mu \pm \sigma\equal {:.1f}\%$".format(stdDev1p * 100))
        plt.grid()
        plt.title(r"Apogee PDF (ID: {}, $\mu = {:.1f}, \sigma = {:.1f}$)".format(self.id, mn, stdDev))
        plt.xlabel("Apogee [m]")
        plt.plot(mn, 0, 'x', color = (1, 0, 0), label = r"$\mu$")
        plt.legend()
        plt.show()

    def europa(self):
        apogees = self.apogees[np.nonzero(self.apogees)]
        mn = np.mean(apogees)
        pdf = gaussian_kde(apogees)
        x = np.linspace(np.min(apogees) - 200, np.max(apogees) + 200, 1000)

        color = (np.random.choice([0, 0.5, 1]), np.random.choice([0, 0.5, 1]), np.random.choice([0, 0.5, 1]))
        while color == (0, 0, 0) or color == (1, 1, 1):
            color = (np.random.choice([0, 0.5, 1]), np.random.choice([0, 0.5, 1]), np.random.choice([0, 0.5, 1]))
        plt.plot(x, pdf(x), color=color, label = self.id + "(mass: {:.2f}, cg: {:.2f})".\
        format(self.rocket.getMass(0), -self.rocket.getCOM(0)[0]))
        plt.plot(mn, 0, 'x', color=color)

    def nigeria(self, x):
        apogees = self.apogees[np.nonzero(self.apogees)]
        pdf = gaussian_kde(apogees)
        cdf = np.zeros(len(x))
        for i in range(len(x)):
            cdf[i] = np.trapz(pdf(x)[0:i], x[0:i])
        return cdf

    @staticmethod
    def africa(x, lower, upper, goTrgt, goVal):
        n = len(lower)
        ids = []
        ps = np.zeros(n, dtype=interp1d)
        for i in range(n):
            ids += [lower[i].id + "/" + upper[i].id]
            psVals = lower[i].nigeria(x) * (1 - upper[i].nigeria(x))
            ps[i] = interp1d(x, psVals)
        goUpper = np.zeros(len(x))
        goUpper.fill(goVal)
        plt.fill_between(x, goUpper, color=(0.8, 0.8, 0.8))
        for i in range(n):
            if ps[i](goTrgt) < goVal:
                color = (1, 0, 0)
                status = "no go"
            else:
                color = (0, 0, 0)
                status = "go"
            plt.plot(x, ps[i](x), label=ids[i] + " ({})".format(status))
            plt.plot(goTrgt, ps[i](goTrgt), 'x', color=color)
        plt.ylabel("probabillity [decimal]")
        plt.xlabel("apogee [m]")
        plt.title("MonteCarlo CDFs (go={}, target={})".format(goVal, goTrgt))
        plt.grid()
        plt.legend()
        plt.show()

    @staticmethod
    def africaBar(x, lower, upper, goTrgt, goVal):
        n = len(lower)
        ids = []
        ps = np.zeros(n, dtype=interp1d)
        for i in range(n):
            ids += [lower[i].id + "/" + upper[i].id]
            psVals = lower[i].nigeria(x) * (1 - upper[i].nigeria(x))
            ps[i] = interp1d(x, psVals)
        goUpper = np.zeros(len(x))
        goUpper.fill(goVal)
        probs = []
        colors = []
        for i in range(n):
            prob = ps[i](goTrgt)
            if prob < goVal:
                color = (1, 0, 0)
                status = "no go"
            else:
                color = (0, 1, 0)
                status = "go"
            probs += [prob]
            colors += [color]
        plt.ylabel("probabillity [decimal]")
        plt.xlabel("Monte Carlo IDs")
        plt.title("Go / No go (go={}, target={})".format(goVal, goTrgt))
        xpos = np.arange(len(ids))
        plt.bar(xpos, probs, color=colors)
        plt.xticks(xpos, ids)
        plt.setp(plt.gcf().gca().get_xticklabels(), rotation=15, horizontalalignment='right')
        plt.show()


    def __str__(self):
        return "id:".ljust(20) + self.id + '\n' + \
        "n:".ljust(20) + str(self.n) + '\n' + \
        "Thrust deviation:".ljust(20) + str(self.thrustDev) + '\n' + \
        "Thrust frequency:".ljust(20) + str(self.thrustFreq) + '\n' + \
        "Drag deviation:".ljust(20) + str(self.dragDev) + 2*'\n' + \
        self.wind.__str__() + '\n'
