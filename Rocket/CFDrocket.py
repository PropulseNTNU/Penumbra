"""
The definition of the CFD rocket

Version: WIP
Last edit: 06.11.2018

--Propulse NTNU--
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
from Rocket import Motor, find_parameter

class CFDrocket:

	def __init__(self, *args):
		print('Initializing rocket..')
		self.__initMass = args[0]
		self.__initInertiaMatrix = args[1]
		self.__initCOM = args[2]
		self.__length = args[3]
		self.__motor = args[4]

		self.__mass = self.__initMass
		self.__inertiaMatrix = self.__initInertiaMatrix
		self.__COM = self.__initCOM


		# Interpolation of forces
		freeAirStreamSpeeds = args[5]
		AoAarray = args[6]
		dragArray = args[7]
		liftArray = args[8]
		momentsArray = args[9]

		print('\tInterpolating the drag force..')
		self.__Dragforce = interp2d(freeAirStreamSpeeds, AoAarray, dragArray)
		print('\tInterpolating the lift force..')
		self.__Liftforce = interp2d(freeAirStreamSpeeds, AoAarray, liftArray)
		print('\tInterpolating the moment about COM..')
		self.__MomentAboutCOM = interp2d(freeAirStreamSpeeds, AoAarray, momentsArray)

		# Center of pressure
		print('\tCalculating the ')

	def getMotor(self):
		return self.__motor

	def getMass(self, t):
		return self.__mass

	def getInertiaMatrix(self, t):
		return self.__inertiaMatrix

	def getCOM(self, t):
		return self.__COM

	# Aero dynamics
	def getDrag(self, speed, AoA):
		return self.__Dragforce(speed, AoA)

	def getLift(self, speed, AoA):
		return self.__Liftforce(speed, AoA)

	def getCOP(self, AoA):
		pass

	# auxiliary
	def plot(self):
		pass

	@staticmethod
	def from_file(initFile, dragFile, liftFile, momentsFile, path_to_file=''):
		"""
		Create an instance of CFDrocket by reading a file with the following
		content:
				- initial mass
				- initial Inertia matrix (about principle axes)
				- initial Center of Mass
				- rocket length
				- rocket motor (a proper file with specs of motor of interest)

				example: initFile.dot

				initMass = value of initial mass
				initMOI = Ixx, Iyy, Izz
				initCOM = x-value of initial COM
				length = value of rocket total length
				 motor = myMotor.dot

		:param initFile: The file with content specified above
		:param dragFile: The CFD file with drag forces
		:param liftFile: The CFD file with lift forces
		:param momentsFile: The CFD file with aero moments about COM.
		:param path_to_file: The path to the 'initFile' relative to current path

		return: A rocket instance with specs from file and CFD.
		"""
		path = path_to_file + initFile
		initMass = find_parameter(path, 'initMass')
		initMOI = find_parameter(path, 'initMOI')
		initMOI = np.diag(float(np.array([x.strip() for x in initMOI.split(',')])))
		initCOM = find_parameter(path, 'initCOM')
		length = find_parameter(path, 'length')
		motor = Motor.from_file(find_parameter(path, 'motor'))

		# aero dynamics


		return CFDrocket(initMass, initMOI, initCOM, length, motor)

def unwrapCFD_files(dragFile, liftFile, momentsFile):
	pass