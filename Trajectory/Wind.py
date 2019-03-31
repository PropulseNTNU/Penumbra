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
        if type(speed0) == list:
            self.speed0 = np.random.uniform(speed0[0], speed0[1])
        else:
            self.speed0 = speed0
        if type(speed0) == list:
            self.direction = np.random.uniform(direction[0], direction[1])
        else:
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
    def __init__(self, alt0, speed0, t, direction):
        self.initParams = [alt0, speed0, t, direction]
        engWind.__init__(self, alt0, speed0, direction)
        self.t = t
        self.freq = 20
        timeDomain = np.linspace(0, t, self.freq * t)
        varience = np.array([np.random.normal(scale = 0.5)\
        for timeDomain in timeDomain])

        self.Uv = interp1d(timeDomain, varience)

    def refresh(self):
        self.__init__(*self.initParams)

    def __str__(self):
        outString = "White wind\n" + "-"*16 +\
        "\nalpha: {:.3f}\nalt0: {:.3f}\nspeed0: {:.3f}\ndirection: {:.3f}"\
        .format(self.alpha, self.alt0, self.speed0, self.direction)
        return outString

    def getMagnitude(self, alt, timeDomain):
        n = self.Uv(timeDomain)
        if alt > 0:
            return self.speed0 * (alt / self.alt0) ** (1 / self.alpha) + n
        else:
            return 0

    def getWindVector(self, alt, timeDomain):
        mag = self.getMagnitude(alt, timeDomain)
        direction = self.direction
        return np.array([mag*np.cos(direction), mag*np.sin(direction), 0])

class pinkWind(whiteWind):
    def __init__(self, alt0, speed0, t, direction = np.random.uniform(0, 2 * np.pi)):
        self.initParams = [alt0, speed0, t, direction]
        alpha = 5/3
        whiteWind.__init__(self, alt0, speed0, t, direction)
        self.t = t
        self.freq = 20
        timeDomain = np.linspace(0, t, self.freq * t)

        n = int(np.ceil(self.freq*t))
        x = np.zeros(n)
        for i in range(2, n):
            a1 = (-alpha/2)
            a2 = (1 - alpha/2)*(a1 / 2)
            x[i] = np.random.normal() - a1 * x[i-1] - a2 * x[i - 2]

        self.Uv = interp1d(timeDomain, x)

    def __str__(self):
        outString = "Pink wind\n" + "-"*16 +\
        "\nalpha: {:.3f}\nalt0: {:.3f}\nspeed0: {:.3f}\ndirection: {:.3f}"\
        .format(self.alpha, self.alt0, self.speed0, self.direction)
        return outString
