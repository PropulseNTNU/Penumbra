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
from mpl_toolkits.mplot3d import Axes3D
from matplotlib import cm

def empty_filter(x):
	if x == '':
		return False
	else:
		return True

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

def step2alphaspeed(T, alpha_max, step, delta_v, v0):
	f = 1 / T
	alpha = (alpha_max/2)*(1+np.sin((T/(T-1))*(step - np.floor(step * f))*np.pi*f-np.pi/2))
	v = v0 + np.floor(step * f) * delta_v
	return alpha, v

def plot_step2alphaspeed(T, alpha_max, delta_v, v0, rng):
	alpha_list = [step2alphaspeed(T, alpha_max, i, delta_v, v0)[0] for i in range(rng)]
	v_list = [step2alphaspeed(T, alpha_max, i, delta_v, v0)[1] for i in range(rng)]

	fig, (ax1, ax2) = plt.subplots(2, 1, sharex = True)
	ax1.plot(v_list, 'o', color='r')
	ax2.plot(alpha_list, 'o', color='b')
	ax1.grid()
	ax2.grid()
	ax1.set_title("speed")
	ax2.set_title("AoA")
	fig.suptitle("step to AoA and speed")
	ax2.set_xlabel("step")
	ax2.set_ylabel("AoA [deg]")
	ax1.set_ylabel("Speed [mach]")
	plt.show()

def make2d(report, c, T, idx):
	arr = report[idx]
	arr = arr[:c]
	arr = np.reshape(arr, (-1, T))
	for i in range(len(arr)):
		if i % 2 == 0:
			arr[i] = np.flip(arr[i])
	return arr

def unwrap_report(file, T, alpha_max, delta_v, v0):
	f = 1/T
	report = np.loadtxt(file, skiprows = 3).transpose()[1:]
	c = len(report.transpose())
	c -= c % T

	drag = make2d(report, c, T, 0)
	lift = make2d(report, c, T, 1)
	moix =  make2d(report, c, T, 2)
	moiy = make2d(report, c, T, 3)
	moiz = make2d(report, c, T, 4)

	# Speed axis seems to be working fine regardless of T
	# Alpha axis acts weird when given other integers than 5
	alpha_axis = alpha_max * np.sin(np.linspace(0, T, T) * np.pi * f/2)
	speed_axis =  v0 + (np.linspace(0, c / T - 1, c // T) * delta_v)
	return drag, lift, mom_x, mom_y, mom_z, alpha_axis, speed_axis

def plot_report(file, T, alpha_max, delta_v, v0):
	drag, lift, moix, moiy, moiz, alpha_axis, speed_axis = unwrap_report(file, T, alpha_max, delta_v, v0)
	fig, (ax1, ax2) = plt.subplots(1, 2, sharey = True)
	fig.suptitle("Force report")
	ax1.imshow(drag, cmap='inferno_r')
	ax2.imshow(lift, cmap='inferno')
	ax1.set_title("Drag")
	ax1.set_ylabel("Speed [unknown]")
	ax1.set_xlabel("AoA [unknown]")
	ax2.set_title("Lift")
	ax2.set_ylabel("Speed [unknown]")
	ax2.set_xlabel("AoA [unknown]")
	print(alpha_axis)
	print(speed_axis)
	plt.show()

plot_report("full-report.out", 5, 10, 0.1, 0.1)
