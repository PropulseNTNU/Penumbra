import numpy as np
from scipy.interpolate import interp1d

class weibull:
    def __init__():
        pass

class nullWind:
    def __init__(self):
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
        self.speed0 = speed0
        self.direction = direction

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

class whiteWind(engWind):
    def __init__(self, alt0, speed0, direction, timeDomain):
        engWind.__init__(self, alt0, speed0, direction)
        self.timeDomain = timeDomain
        self.freq = 20
        t = np.linspace(0, timeDomain, self.freq * timeDomain)
        varience = np.array([np.random.normal(scale = 0.5)\
        for t in t])

        self.Uv = interp1d(t, varience)


    def __str__(self):
        outString = "White wind\n" + "-"*16 +\
        "\nalpha: {:.3f}\nalt0: {:.3f}\nspeed0: {:.3f}\ndirection: {:.3f}"\
        .format(self.alpha, self.alt0, self.speed0, self.direction)
        return outString

    def getMagnitude(self, alt, t):
        n = self.Uv(t)
        if alt > 0:
            return self.speed0 * (alt / self.alt0) ** (1 / self.alpha) + n
        else:
            return 0

    def getWindVector(self, alt, t):
        mag = self.getMagnitude(alt, t)
        direction = self.direction
        return np.array([mag*np.cos(direction), mag*np.sin(direction), 0])
