import numpy as np

class weibull:
    def __init__():
        pass

class nullWind:
    def __init__(self):
        pass

    def getWindVector(self, alt):
        return np.zeros(3)

    def __str__(self):
        outString = "nullWind\n"  + "-"*16 + "\nMy purpouse is to act as a\
default when no windObj is specified"
        return outString

class engWind():
    def __init__(self, alt0, speed0, direction):
        self.__alpha = 1 / 0.143
        self.__alt0 = alt0
        self.__speed0 = speed0
        self.__direction = direction

    def __str__(self):
        outString = "Engineering wind\n" + "-"*16 +\
        "\nalpha: {:.3f}\nalt0: {:.3f}\nspeed0: {:.3f}\ndirection: {:.3f}"\
        .format(self.__alpha, self.__alt0, self.__speed0, self.__direction)
        return outString

    def getMagnitude(self, alt):
        if alt > 0:
            return self.__speed0 * (alt / self.__alt0) ** (1 / self.__alpha)
        else:
            return 0

    def getWindVector(self, alt):
        mag = self.getMagnitude(alt)
        direction = self.__direction
        return np.array([mag*np.cos(direction), mag*np.sin(direction), 0])
