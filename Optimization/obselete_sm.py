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
            apogee = np.max(-Trajectory.calculateTrajectoryWithAirbrakes(self.rocket,\
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

    def stochasticShoot(self, t0, dt, targetApogee, tol):
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
