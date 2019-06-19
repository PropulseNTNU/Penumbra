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
# Import of external libs
import numpy as np
import matplotlib.pyplot as plt
# Import of internal libs
import Wind
from Rocket1 import RocketSimple
import TrajectoryWithBrakes as Trajectory
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

class ShootingMethod:
    def __init__ (self, t, dt, target, mc, tol):
        self.t = t
        self.dt = dt
        self.target = target
        self.mc = mc
        self.tol = tol
        self.apogee = 0
        self.mc = mc
        self.mc.setVerbose(False)

    def run(self):
        while (np.abs(self.apogee - self.target) > self.tol):
            i = np.sign(self.apogee - self.target)

            self.mc.flush()
            self.mc.setTbrakes(self.t)
            self.mc.run()

            self.showProgess()

            self.apogee = self.mc.getMeanApogee()
            j = np.sign(self.apogee - self.target)

            if i != j: self.dt /= 2

            if self.apogee < self.target - self.tol: self.t += self.dt
            if self.apogee > self.target + self.tol: self.t -= self.dt

        print(self.genTrajectory())

    def genTrajectory(self):
        traj = self.mc.getTrajectory(self.target, self.tol)
        pos = -traj[1].T[2]
        vel = -traj[4].T[2]
        ap_i = np.where(vel < 0)[0][0]
        apogee = np.max(pos)
        pos = pos[0:ap_i]
        vel = vel[0:ap_i]
        vel_func = interp1d(pos, vel)
        pos_of_interest = np.arange(10, int(apogee), 1)
        vel_of_interest = vel_func(pos_of_interest)

        out = open("out.txt", 'w')
        for i in range(len(pos_of_interest)):
            out.write(str(pos_of_interest[i]) + '\t' + str(vel_of_interest[i]) + '\n')
        out.close()
        return apogee

    def showProgess(self):
        print(self.t, self.dt, self.apogee)


class MonteCarlo:
    def __init__(self, n, rocket, params, thrustFreq = 0, thrustDev = 0,\
    dragDev = 0, wind = Wind.nullWind(), doSave = True, verbose = True,\
    Tbrakes = 1e30, Cbrakes_in = 0, target = 3048, type = 'default'):
        self.n = n
        self.type = type
        self.target = target
        self.deltas = np.zeros(n)
        self.params = params
        self.rocket = rocket
        self.rocket.getMotor().setStochasticParams(thrustFreq, thrustDev)
        self.rocket.setDragDev(dragDev)
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
        self.trajectories = [0]*n
        if type == 'default':
            self.apogees = np.zeros(n)
            self.run = self.run_default
        if type == 'gng':
            self.abmax_apogees = np.zeros(n)
            self.abmin_apogees = np.zeros(n)
            self.Tbrakes = self.rocket.getMotor().getBurnTime()
            self.run = self.run_gng
            self.successArr = np.zeros(n)
            self.success = 0
        if self.verbose: print("Monte")

    def showProgess_default(self):
        scaledown = 1.48
        fullbar = int(100 // scaledown)
        progress = int((self.i)/self.n * 100)
        pbar = int(progress // scaledown)
        meanApogee = np.mean(self.apogees[np.nonzero(self.apogees)])
        meanDelta = np.mean(self.deltas[np.nonzero(self.deltas)])
        timeLeft = meanDelta * (self.n - self.i)
        stdDevApogee = np.std(self.apogees[np.nonzero(self.apogees)])

        os.system("cls")
        print(\
        "Iterations".ljust(20) + \
        "Progress".ljust(20) + \
        "Time left".ljust(20) + \
        "Mean apogee".ljust(20) + \
        "Std. deviation".ljust(20))

        print('-'*71)

        print(\
        "{} / {}".format(self.i, self.n).ljust(20) + \
        "{}%".format(progress).ljust(20) + \
        "{}".format(str(dt.timedelta(seconds=int(timeLeft)))).ljust(20) + \
        "{}m".format(round(meanApogee, 2)).ljust(20) + \
        "{}m".format(round(stdDevApogee, 2)).ljust(20))
        print('\n' + "  " + pbar*chr(9619) + (fullbar - pbar)*chr(9617) + '\n')

        print(self)

    def showProgess_gng(self):
        scaledown = 1.48
        fullbar = int(100 // scaledown)
        progress = int((self.i)/self.n * 100)
        pbar = int(progress // scaledown)
        abmax_meanApogee = np.mean(self.abmax_apogees[np.nonzero(self.abmax_apogees)])
        abmin_meanApogee = np.mean(self.abmin_apogees[np.nonzero(self.abmin_apogees)])
        meanRange = np.abs(abmin_meanApogee - abmax_meanApogee)
        meanDelta = np.mean(self.deltas[np.nonzero(self.deltas)])
        timeLeft = meanDelta * (self.n - self.i)
        success = self.success

        os.system("cls")
        print(\
        "Iterations".ljust(20) + \
        "Progress".ljust(20) + \
        "Time left".ljust(20) + \
        "Mean apogee (ab = 0)".ljust(20) + \
        "Mean apogee (ab = C)".ljust(20) + \
        "Mean range".ljust(20) + \
        "Success".ljust(20))

        print('-'*141)

        print(\
        "{} / {}".format(self.i, self.n).ljust(20) + \
        "{}%".format(progress).ljust(20) + \
        "{}".format(str(dt.timedelta(seconds=int(timeLeft)))).ljust(20) + \
        "{}m".format(round(abmin_meanApogee, 2)).ljust(20) + \
        "{}m".format(round(abmax_meanApogee, 2)).ljust(20) + \
        #"{}m".format(round(self.abmin_apogees[self.i - 1], 2)).ljust(20) + \
        #"{}m".format(round(self.abmax_apogees[self.i - 1], 2)).ljust(20) + \
        "{}m".format(round(meanRange, 2)).ljust(20) + \
        "{}%".format(round(success*100, 2)).ljust(20))
        print('\n' + "  " + pbar*chr(9619) + (fullbar - pbar)*chr(9617) + '\n')

        print(self)

    def print_gng(self):
        print(str(round(self.meanApogee, 2)) + "+-" + str(round(self.meanRange/2, 2)))

    def run_default(self):
        global interrupted
        #os.system("cls")
        def signal_handler(signal, frame):
            global interrupted
            print("Process safely interrupted")
            interrupted = True
            return (self.i)/self.n * 100

        signal.signal(signal.SIGINT, signal_handler)

        while self.i < self.n and not interrupted:
            T0 = dt.datetime.now()
            self.rocket.refresh()
            self.wind.refresh()

            trajectory = Trajectory.calculateTrajectoryWithAirbrakes(self.rocket,\
            *self.params[0:4], windObj = self.wind,\
            dragDeviation = self.dragDev, Tbrakes = self.Tbrakes,\
            Cbrakes_in = self.Cbrakes_in)


            self.apogees[self.i] = np.max(-trajectory[1][:,2])
            self.trajectories[self.i] = trajectory

            if self.doSave:
                self.save()

            Tf = dt.datetime.now()
            self.deltas[self.i] = (Tf - T0).total_seconds()

            self.i += 1

            if self.verbose: self.showProgess_default()

        interrupted = False
        if self.verbose: print("Process terminated, id: {}".format(self.id))
        return (self.i)/self.n * 100

    def run_gng(self):
        global interrupted
        #os.system("cls")
        def signal_handler(signal, frame):
            global interrupted
            print("Process safely interrupted")
            interrupted = True
            return (self.i)/self.n * 100

        signal.signal(signal.SIGINT, signal_handler)

        while self.i < self.n and not interrupted:
            T0 = dt.datetime.now()
            self.rocket.refresh()
            self.wind.refresh()

            abmin_trajectory = Trajectory.calculateTrajectoryWithAirbrakes(self.rocket,\
            *self.params[0:4], windObj = self.wind,\
            dragDeviation = self.dragDev, Tbrakes = 10000,\
            Cbrakes_in = 0)

            abmax_trajectory = Trajectory.calculateTrajectoryWithAirbrakes(self.rocket,\
            *self.params[0:4], windObj = self.wind,\
            dragDeviation = self.dragDev, Tbrakes = self.Tbrakes,\
            Cbrakes_in = self.Cbrakes_in)

            self.abmin_apogees[self.i] = np.max(-abmin_trajectory[1][:,2])
            self.abmax_apogees[self.i] = np.max(-abmax_trajectory[1][:,2])

            if self.abmin_apogees[self.i] > self.target and\
            self.abmax_apogees[self.i] < self.target:
                self.successArr[self.i] = 1

            self.success = np.sum(self.successArr) / (self.i + 1)

            if self.doSave:
                self.save()

            Tf = dt.datetime.now()
            self.deltas[self.i] = (Tf - T0).total_seconds()

            self.i += 1

            if self.verbose: self.showProgess_gng()

        interrupted = False
        if self.verbose: print("Process terminated, id: {}".format(self.id))
        self.meanRange = np.abs(np.mean(self.abmin_apogees - self.abmax_apogees))
        self.abmax_meanApogee = np.mean(self.abmax_apogees)
        self.abmin_meanApogee = np.mean(self.abmin_apogees)
        self.meanApogee = np.mean([self.abmax_meanApogee, self.abmin_meanApogee])
        return (self.i)/self.n * 100

    def getMeanDelta():
        return self.meanDelta

    def getDeltas():
        return self.deltas

    def getMeanApogee(self):
        return np.mean(self.apogees[0:self.i])

    def setVerbose(self, bol):
        self.verbose = bol

    def setTbrakes(self, t):
        self.Tbrakes = t

    def getTrajectory(self, target, tol):
        for traj in self.trajectories[0:self.i]:
            if np.abs(np.max(-traj[1][:,2]) - target) < tol: return traj

    def genTrajectory(self):
        traj = self.getTrajectory(self.getMeanApogee(), 1)
        pos = -traj[1].T[2]
        vel = -traj[4].T[2]
        ap_i = np.where(vel < 0)[0][0]
        apogee = np.max(pos)
        pos = pos[0:ap_i]
        vel = vel[0:ap_i]
        vel_func = interp1d(pos, vel)
        pos_of_interest = np.arange(10, int(apogee), 1)
        vel_of_interest = vel_func(pos_of_interest)

        out = open("out.txt", 'w')
        for i in range(len(pos_of_interest)):
            out.write(str(pos_of_interest[i]) + '\t' + str(vel_of_interest[i]) + '\n')
        out.close()
        return apogee

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

    def europa(self, c=0):
        apogees = self.apogees[np.nonzero(self.apogees)]
        mn = np.mean(apogees)
        pdf = gaussian_kde(apogees)
        x = np.linspace(np.min(apogees) - 200, np.max(apogees) + 200, 1000)

        if c == 0:
            color = (np.random.choice([0, 0.5, 1]), np.random.choice([0, 0.5, 1]), np.random.choice([0, 0.5, 1]))
            while color == (0, 0, 0) or color == (1, 1, 1):
                color = (np.random.choice([0, 0.5, 1]), np.random.choice([0, 0.5, 1]), np.random.choice([0, 0.5, 1]))
        else:
            color = c
        plt.plot(x, pdf(x), color=color, label = self.id)
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
        "Drag deviation:".ljust(20) + str(self.dragDev) + '\n' + \
        "Type:".ljust(20) + str(self.type) + 2*'\n' + \
        self.wind.__str__() + '\n'
