"""
The definition of the CFD rocket

Version: WIP
Last edit: 11.11.2018

--Propulse NTNU--
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
from Rocket.Rocket1 import Motor
from Rocket.lib.File_utilities import find_parameter, unwrap_report

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
		dragMatrix = args[7]
		liftMatrix = args[8]
		momentsMatrix_x = args[9]
		momentsMatrix_y = args[10]
		momentsMatrix_z = args[11]

		print('\tInterpolating the drag force..')
		self.__Dragforce = interp2d(freeAirStreamSpeeds, AoAarray, dragMatrix)
		print('\tInterpolating the lift force..')
		self.__Liftforce = interp2d(freeAirStreamSpeeds, AoAarray, liftMatrix)
		print('\tInterpolating the moment about COM (x-component)..')
		self.__MomentAboutCOM_X = interp2d(freeAirStreamSpeeds, AoAarray, momentsMatrix_x)
		print('\tInterpolating the moment about COM (y-component)..')
		self.__MomentAboutCOM_Y = interp2d(freeAirStreamSpeeds, AoAarray, momentsMatrix_y)
		print('\tInterpolating the moment about COM (z-component)..')
		self.__MomentAboutCOM_Z = interp2d(freeAirStreamSpeeds, AoAarray, momentsMatrix_z)
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
		return self.__length

	# Aero dynamics
	def getDrag(self, speed, AoA):
		"""
		:param speed: [float] the air speed relative to rocket [m/s]
		:param AoA: [float] the angle of attack [rad]

		:return: [float] The Drag on rocket attacking in COP [N]
		"""
		return self.__Dragforce(speed, AoA)

	def getLift(self, speed, AoA):
		"""
		:param speed: [float] the air speed relative to rocket [m/s]
		:param AoA: [float] the angle of attack [rad]

		:return: [float] The Lift on rocket attacking in COP [N]
		"""
		return self.__Liftforce(speed, AoA)

	def getMomentAboutCOM_X(self, speed, AoA):
		"""
		:param speed: [float] the air speed relative to rocket [m/s]
		:param AoA: [float] the angle of attack [rad]

		:return: [float] The total moment on rocket about COM (x-component) [Nm]
		"""
		return self.__MomentAboutCOM_X(speed, AoA)

	def getMomentAboutCOM_Y(self, speed, AoA):
		"""
		:param speed: [float] the air speed relative to rocket [m/s]
		:param AoA: [float] the angle of attack [rad]

		:return: [float] The total moment on rocket about COM (y-component) [Nm]
		"""
		return self.__MomentAboutCOM_Y(speed, AoA)

	def getMomentAboutCOM_Z(self, speed, AoA):
		"""
		:param speed: [float] the air speed relative to rocket [m/s]
		:param AoA: [float] the angle of attack [rad]

		:return: [float] The total moment on rocket about COM (z-component) [Nm]
		"""
		return self.__MomentAboutCOM_Z(speed, AoA)


	def getCOP(self, speed, AoA):
		"""
		:param speed: [float] the air speed relative to rocket [m/s]
		:param AoA: [float] the angle of attack [rad]

		:return: [float] The position of COP relative to nose tip [m]
		"""
		Fd = self.getDrag(speed, AoA)
		Fl= self.getLift(speed, AoA)
		#M = self.getMomentAboutCOM(speed, AoA)
		COM_0 = self.getCOM(0)

		return COM_0 - M/(Fl*np.cos(AoA) + Fd*np.sin(AoA))

	def getInertiaMatrix(self, t):
		"""
		:param t: [float] point in time [s]
		:return [np.array, 3x3 matrix] The inertia matrix at time t [kgm^2]

		NOTE: Currently assuming that fuel is burning radially.
		"""
		#TODO Add inertia for axial burning fuel.
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
		return I0 + rInitMass*np.diag([0, deltaR**2, deltaR**2]) + (mInertia - mInitInertia) + (mMass - mInitMass)*np.diag([0, deltaM**2, deltaM**2])

	# auxiliary
	def plot(self):
		#Plots of motor performance
		self.__motor.plotPerformance()

		#Plots of rocket characteristics (during burn phase)
		time = self.__time
		burnTime = self.__motor.getBurnTime()
		print('Creating plot of COM and mass...')
		COM = np.array([self.getCOM(t) for t in time])
		mass = np.array([self.getMass(t) for t in time])
		ax1 = plt.subplot(211, xlabel='time [s]', ylabel='[cm]')
		ax1.plot(time, 100*COM, label='COMx(t)', lw=2, c='r')
		ax1.set_title('COM of rocket during burn phase of %1.1f s' % burnTime)
		ax1.grid()

		ax2 = plt.subplot(212, xlabel='time [s]', ylabel='[kg]')
		ax2.plot(time, mass, label='mass(t)', lw=2, c='b')
		ax2.set_title('Mass of rocket during burn phase of %1.1f s' % burnTime)
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
		ax1.plot(time, Ixx*1e4, label='Ixx(t)', lw=2, c='r') #Multiply by 10,000 to get correct unit
		ax1.set_title('Ixx of rocket during burn phase of %1.1f s' % burnTime)
		ax1.grid()

		ax2 = plt.subplot(212, xlabel='time [s]', ylabel='[kgcm²]')
		ax2.plot(time, Iyy*1e4, label='I(t)', lw=2, c='b')  #Multiply by 10,000 to get correct unit
		ax2.set_title('Iyy and Izz during burn phase of %1.1f s' % burnTime)
		ax2.grid()
		plt.subplots_adjust(hspace=0.5)
		plt.show()
		print('Plotting done!\n')

	@staticmethod
	def from_file(initFile, sampleReport, path_to_file=''):
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
		initMass = find_parameter(path, 'initial_mass')
		initMOI = find_parameter(path, 'initial_moi')
		initMOI = np.diag(np.array([x.strip() for x in initMOI.split(',')]).astype(float))
		initCOM = find_parameter(path, 'initial_com')
		length = find_parameter(path, 'length')
		motor = Motor.from_file(path_to_file + find_parameter(path, 'motor'))

		#TODO Implement function to read CFD files.
		path = path_to_file + sampleReport
		T = find_parameter(path, 'period')
		alpha_max = find_parameter(path, 'alpha_max')
		delta_v = find_parameter(path, 'delta_v')
		v0 = find_parameter(path, 'v0')
		drag, lift, mom_x, mom_y, mom_z, alpha, air_speed = unwrap_report(path,
															int(T), float(alpha_max), float(delta_v), float(v0))

		return Rocket(float(initMass), initMOI, float(initCOM), float(length), motor, air_speed, alpha, drag, lift, mom_x, mom_y, mom_z)
