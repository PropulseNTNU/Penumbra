import numpy as np
from scipy.interpolate import interp1d

class nullWind:
    def __init__(self):
        pass

    def refresh(self):
        pass

    def getWindVector(self, alt, *args):
        return np.zeros(3)

    def __str__(self):
        outString = "nullWind\n"  + "-"*16 + "\nMy purpouse is to act as a\
default when no windObj is specified"
        return outString

class engWind:
    def __init__(self, alt0, speed0, direction):
        self.alpha = 1 / 0.143
        self.alt0 = alt0
        if type(speed0) == list:
            self.speed0 = np.random.uniform(speed0[0], speed0[1])
            self.speed0Array = speed0
        else:
            self.speed0 = speed0
            self.speed0Array = [speed0]

        if type(direction) == list:
            self.direction = np.random.uniform(direction[0], direction[1])
            self.directionArray = direction
        else:
            self.direction = direction
            self.directionArray = [direction]

    def __str__(self):
        outString = "Engineering wind\n" + "-"*16 +\
        "\nalpha: {:.3f}\nalt0: {:.3f}\nspeed0: {:.3f}\ndirection: {:.3f}"\
        .format(self.alpha, self.alt0, self.speed0, self.direction)
        return outString

    def getMagnitude(self, alt):
        if alt > 0:
            return self.speed0 * (alt / self.alt0) ** (1 / self.alpha)
        else:
            return 0

    def getWindVector(self, alt, *args):
        mag = self.getMagnitude(alt)
        direction = self.direction
        return np.array([mag*np.cos(direction), mag*np.sin(direction), 0])

class pinkWind(engWind):
    def __init__(self, t, speed0, alt0 = 10, intensity = 0.1, direction = 0):
        self.initParams = [t, speed0, alt0, intensity, direction]
        engWind.__init__(self, alt0, speed0, direction)
        # To avoid extrapolation errors
        t += 1
        alpha = 5/3
        self.t = t
        self.freq = 20
        timeDomain = np.linspace(0, t, self.freq * t)

        if type(intensity) == list:
            self.I = np.random.uniform(intensity[0], intensity[1])
            self.IArray = intensity
        else:
            self.I = intensity
            self.IArray = [intensity]

        n = int(np.ceil(self.freq*t))
        x = np.zeros(n)
        for i in range(2, n):
            a1 = (-alpha/2)
            a2 = (1 - alpha/2)*(a1 / 2)
            x[i] = np.random.normal() - a1 * x[i-1] - a2 * x[i - 2]

        self.Uv = interp1d(timeDomain, x)

    def refresh(self):
        self.__init__(*self.initParams)

    def getMagnitude(self, alt, timeDomain):
        if alt > 0:
            U = self.speed0 * (alt / self.alt0) ** (1 / self.alpha)
            n = self.Uv(timeDomain) * (self.I * U)
            return U + n
        else:
            return 0

    def getWindVector(self, alt, timeDomain):
        mag = self.getMagnitude(alt, timeDomain)
        direction = self.direction
        return np.array([mag*np.cos(direction), mag*np.sin(direction), 0])

    def __str__(self):
        outString = "Pink wind\n" + "-"*16 +\
        "\nalpha: {:.3f}\nalt0: {:.3f}\nspeed0: {}\ndirection: {}\nintensity: {}"\
        .format(self.alpha, self.alt0, self.speed0Array, self.directionArray, self.IArray)
        return outString
