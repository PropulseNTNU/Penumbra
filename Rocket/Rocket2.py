"""
The definition of the CFD rocket

Version: WIP
Last edit: 16.11.2018

--Propulse NTNU--
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
from Rocket.Rocket1 import Motor
from Rocket.lib.File_utilities import find_parameter, unwrap_report1, unwrap_report2

# Define some things for plotting
font = {'family': 'sans-serif', 'weight': 'bold', 'size': 16}
plt.rc('font', **font)
plt.rcParams['text.latex.preamble'] = [r'\boldmath']


class Rocket:

	def __init__(self, *args):
		print('Initializing rocket..')
		self.__initMass = args[0]
		self.__initInertiaMatrix = args[1]
		self.__initCOM = args[2]
		self.__length = args[3]
		self.__motor = args[4]

		# for plotting, array of time during burn phase

		self.__time = np.arange(0, self.__motor.getBurnTime(), 5e-4)

		# Interpolation of forces and moments
		freeAirStreamSpeeds = args[5]
		AoAarray = args[6]
		dragArray = args[7]
		liftArray = args[8]
		momentsArray_CG = args[9]

		print('\tInterpolating the drag force..')
		self.__Dragforce = interp2d(AoAarray, freeAirStreamSpeeds, dragArray)
		print('\tInterpolating the lift force..')
		self.__Liftforce = interp2d(AoAarray, freeAirStreamSpeeds, liftArray)
		print('\tInterpolating the moment about COM (component normal to aerodynamic plane)..')
		self.__MomentAboutCOM = interp2d(AoAarray, freeAirStreamSpeeds, momentsArray_CG)

		# Done
		print('Rocket initialized!\n')

	def getMotor(self):
		return self.__motor

	def getMass(self, t):
		return self.__initMass + self.__motor.getMass(t) - self.__motor.getMass(0)

	def getCOM(self, t):
		"""
		:param t: [float] at time t [sec]

		:return: [float] the position of COM relative to nose tip [m]
		"""
		mass = self.getMass(t)
		initMass = self.__initMass
		initCOM = self.__initCOM
		rocketLength = self.__length
		motorHeight = self.getMotor().getLength()
		initMotorMass = self.getMotor().getMass(0)
		motorMass = self.getMotor().getMass(t)

		return 1/mass*(initMass*initCOM + (rocketLength - motorHeight/2)*(initMotorMass - motorMass))

	def getLength(self):
		"""
		:return: the length of the rocket [m]
		"""
		return self.__length

	# Aero dynamics
	def getDrag(self, speed, AoA):
		"""
		:param speed: [float] the air speed relative to rocket [m/s]
		:param AoA: [float] the angle of attack [rad]

		:return: [float] The Drag on rocket attacking in COP [N]
		"""
		return self.__Dragforce(AoA, speed)

	def getLift(self, speed, AoA):
		"""
		:param speed: [float] the air speed relative to rocket [m/s]
		:param AoA: [float] the angle of attack [rad]

		:return: [float] The Lift on rocket attacking in COP [N]
		"""
		return self.__Liftforce(AoA, speed)

	def getMomentAboutCOM(self, speed, AoA):
		"""
		:param speed: [float] the air speed relative to rocket [m/s]
		:param AoA: [float] the angle of attack [rad]

		:return: [float] The total moment on rocket about COM (component normal to aerodynamic plane) [Nm]
		"""
		return self.__MomentAboutCOM(AoA, speed)

	def getCOP(self, speed, AoA):
		"""
		:param speed: [float] the air speed relative to rocket [m/s]
		:param AoA: [float] the angle of attack [rad]

		:return: [float] The position of COP relative to nose tip [m]
		"""
		Fd = self.getDrag(AoA, speed)
		Fl = self.getLift(AoA, speed)
		M = self.getMomentAboutCOM(AoA, speed)
		COM_0 = self.getCOM(0)

		return COM_0 - M/(Fl*np.cos(AoA) + Fd*np.sin(AoA))

	def getInertiaMatrix(self, t):
		"""
		:param t: [float] point in time [s]
		:return [np.array, 3x3 matrix] The inertia matrix at time t [kgm^2]

		NOTE: Currently assuming that fuel is burning radially.
		"""
		# TODO Add inertia for axial burning fuel.
		I0 = self.__initInertiaMatrix
		rInitMass = self.__initMass
		rInitCOM = self.__initCOM
		rCOM = self.getCOM(t)
		mInitMass = self.__motor.getMass(0)
		mInitCOM = self.__motor.getCOM(0)
		mInitInertia = self.__motor.getInertiaMatrix(0)
		mMass = self.__motor.getMass(t)
		mInertia = self.__motor.getInertiaMatrix(t)
		deltaR = rCOM - rInitCOM
		deltaM = rCOM + mInitCOM + self.getLength()
		return I0 + rInitMass*np.diag([0, deltaR**2, deltaR**2]) + (mInertia - mInitInertia) + (
					mMass - mInitMass)*np.diag([0, deltaM**2, deltaM**2])

	# auxiliary
	def plot(self):
		# Plots of motor performance
		self.__motor.plotPerformance()

		# Plots of rocket characteristics (during burn phase)
		time = self.__time
		burnTime = self.__motor.getBurnTime()
		print('Creating plot of COM and mass...')
		COM = np.array([self.getCOM(t) for t in time])
		mass = np.array([self.getMass(t) for t in time])
		ax1 = plt.subplot(211, xlabel='time [s]', ylabel='[cm]')
		ax1.plot(time, 100*COM, label='COMx(t)', lw=2, c='r')
		ax1.set_title('COM of rocket during burn phase of %1.1f s'%burnTime)
		ax1.grid()

		ax2 = plt.subplot(212, xlabel='time [s]', ylabel='[kg]')
		ax2.plot(time, mass, label='mass(t)', lw=2, c='b')
		ax2.set_title('Mass of rocket during burn phase of %1.1f s'%burnTime)
		ax2.grid()
		plt.subplots_adjust(hspace=0.5)
		plt.show()
		# Moment of inertia
		print('Creating plot of Ixx, Iyy and Izz...')
		steps = len(time)
		Ixx = np.zeros(steps)
		Iyy = np.zeros(steps)
		for i in range(steps):
			MOI = np.diagonal(self.getInertiaMatrix(time[i]))
			Ixx[i] = MOI[0]
			Iyy[i] = MOI[1]
		# Izz = Iyy, don't need to create another list for this
		ax1 = plt.subplot(211, xlabel='time [s]', ylabel='[kgcm²]')
		ax1.plot(time, Ixx*1e4, label='Ixx(t)', lw=2, c='r')  # Multiply by 10,000 to get correct unit
		ax1.set_title('Ixx of rocket during burn phase of %1.1f s'%burnTime)
		ax1.grid()

		ax2 = plt.subplot(212, xlabel='time [s]', ylabel='[kgcm²]')
		ax2.plot(time, Iyy*1e4, label='I(t)', lw=2, c='b')  # Multiply by 10,000 to get correct unit
		ax2.set_title('Iyy and Izz during burn phase of %1.1f s'%burnTime)
		ax2.grid()
		plt.subplots_adjust(hspace=0.5)
		plt.show()
		print('Plotting done!\n')

	@staticmethod
	def from_file_without_AoAspeed(initFile, sampleReport, path_to_file=''):
		"""
		Create an instance of CFDrocket by reading some files

		example: initFile.dot as initFile, myMotor.dot as motor and CFD files as drag.dot, lift.dot, moments.dot
							-stored in folder 'myRocket2' (see 'Tests' folder)

				in 'initFile.dot':

				initial_mass = value of initial mass
				initial_moi = Ixx, Iyy, Izz
				initial_com = x-value of initial COM (this one is negative!)
				length = value of rocket total length
				 motor = myMotor.dot

				myRocket = Rocket.from_file('initFile.dot', 'sample-report.dot',
											'myRocket2/')

		:param initFile: The file with content specified above
		:param sampleReport: The CFD file with aero moments about COM.
							Add sampling period, alpha_max, delta_v, v0 values to bottom of sampleReport:
											period = ..
											alpha_max = ..
											.
											.
											.

		:param path_to_file: The path to the files above relative to current path (none by default)

		return: A rocket instance with specs from initFile and CFD.
		"""
		path = path_to_file + initFile
		initMass = find_parameter(path, 'initial_mass')  # in grams
		initMOI = find_parameter(path, 'initial_moi')
		initMOI = np.diag(np.array([x.strip() for x in initMOI.split(',')]).astype(float))  # in g*mm^2
		initCOM = find_parameter(path, 'initial_com') # in millimeters
		length = find_parameter(path, 'length') # in millimeters
		motor = Motor.from_file(path_to_file + find_parameter(path, 'motor'))

		path = path_to_file + sampleReport
		T = find_parameter(path, 'period')
		alpha_max = find_parameter(path, 'alpha_max') # In degrees
		delta_v = find_parameter(path, 'delta_v')  # In mach
		v0 = find_parameter(path, 'v0') # In mach
		alpha, air_speed, drag, lift, moment = unwrap_report1(path,
															  int(T), float(alpha_max), float(delta_v), float(v0))

		return Rocket(float(initMass)/1e3, initMOI/1e9, float(initCOM)/1e3, float(length)/1e3, motor, air_speed,
					  alpha, drag, lift,
					  moment)

	@staticmethod
	def from_file_with_AoAspeed(initFile, sampleReport, path_to_file=''):
		"""
				Create an instance of CFDrocket by reading some files

				example: initFile.dot as initFile, myMotor.dot as motor and CFD files as drag.dot, lift.dot, moments.dot
									-stored in folder 'myRocket2' (see 'Tests' folder)

						in 'initFile.dot':

						initial_mass = value of initial mass
						initial_moi = Ixx, Iyy, Izz
						initial_com = x-value of initial COM (this one is negative!)
						length = value of rocket total length
						 motor = myMotor.dot

						myRocket = Rocket.from_file('initFile.dot', 'sample-report.dot',
													'myRocket2/')

				:param initFile: The file with content specified above
				:param sampleReport: The CFD file with aero moments about COM.
									Add sampling period, alpha_max, delta_v, v0 values to bottom of sampleReport:
													period = ..
													alpha_max = ..
													.
													.
													.

				:param path_to_file: The path to the files above relative to current path (none by default)

				return: A rocket instance with specs from initFile and CFD.
				"""
		path = path_to_file + initFile
		# CAD files are using the units mentioned below
		initMass = find_parameter(path, 'initial_mass')  # in grams
		initMOI = find_parameter(path, 'initial_moi')
		initMOI = np.diag(np.array([x.strip() for x in initMOI.split(',')]).astype(float))  # in g*mm^2
		initCOM = find_parameter(path, 'initial_com')  # in millimeters
		length = find_parameter(path, 'length')  # in millimeters
		motor = Motor.from_file(path_to_file + find_parameter(path, 'motor'))

		path = path_to_file + sampleReport
		alpha, air_speed, drag, lift, moment = unwrap_report2(path)

		return Rocket(float(initMass)/1e3, initMOI/1e9, float(initCOM)/1e3, float(length)/1e3, motor, air_speed,
					  alpha, drag, lift,
					  moment)
