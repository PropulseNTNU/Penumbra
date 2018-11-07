"""
The definition of the ANSYS rocket

Version: WIP
Last edit: 06.11.2018

--Propulse NTNU--
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from Rocket import Motor

class CFDrocket:

	def __init__(self, *args):
		self.__initMass = args[0]
		self.__initInertiaMatrix = np.diag(args[1])
		self.__initCOM = args[2]
		self.__length = args[3]
		self.__motor = args[4]

		self.__mass = self.__initMass
		self.__inertiaMatrix = self.__initInertiaMatrix
		self.__COM = self.__initCOM

	def getMass(self, t):
		return self.__mass

	def getInertiaMatrix(self, t):
		return self.__inertiaMatrix

	def getCOM(self, t):
		return self.__COM

	@staticmethod
	def from_file(initFile):
		pass