"""
Rocket module - The definition of the rocket and its constituent parts

Version: WIP
Last edit: 26.10.2018

--Propulse NTNU--
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import simps
from scipy.interpolate import interp1d

# Geometry types
noseTypes = ['conic', 'dome']


class Nose:

	# Constructor for conic nose
	def __init__(self, noseType, *args):
		# Member variables
		self.__noseType = noseType

		if noseType == noseTypes[0]:  # Conic
			self.__diameter = args[0]
			self.__height = args[1]
			self.__thickness = args[2]
			self.__density = args[3]
		elif noseType == noseTypes[1]:  # Dome
			self.__diameter = args[0]
			self.__thickness = args[1]
			self.__density = args[2]
		print("Nose initialized!")

	def __str__(self):

		if self.__noseType == noseTypes[0]:  # Conic
			D = str(self.__diameter)
			h = str(self.__height)
			d = str(self.__thickness*1e3)
			m = str(round(self.getMass(), 2))
			rho = str(self.__density)
			return "Diameter: " + D + " m\n" + "Height: " + h + " m\n" + "Thickness: " + d + " mm\n" + "Mass: " + m + \
				   " kg\n" + "Density: " + rho + " kgm^-3"

		elif self.__noseType == noseTypes[1]:  # Dome
			D = str(self.__diameter)
			d = str(self.__thickness*1e3)
			m = str(round(self.getMass(), 2))
			rho = str(self.__density)
			return "Diameter: " + D + " m\n" + "Thickness: " + d + " mm\n\n" + "Mass: " + m + \
				   " kg\n" + "Density: " + rho + " kgm^-3"

	# Member functions
	def getVolume(self):
		if self.__noseType == noseTypes[0]:  # Conic
			d = self.__thickness
			R2 = self.__diameter/2  # Outer radius of cone
			H2 = self.__height
			H1 = H2 - d
			R1 = R2*H1/H2
			return np.pi/3*(H2*R2**2 - H1*R1**2)
		elif self.__noseType == noseTypes[1]:  # Dome
			d = self.__thickness
			R2 = self.__diameter/2
			R1 = R2 - d
			return 2/3*np.pi*(R2**3 - R1**3)

	def getCavityVolum(self):
		d = self.__thickness
		R2 = self.__diameter/2  # Outer radius of cone
		H2 = self.__height
		H1 = H2 - d
		R1 = R2*H1/H2
		return np.pi/3*H1*R1**2

	def getLength(self):
		if self.__noseType == noseTypes[0]:
			return self.__height
		elif self.__noseType == noseTypes[1]:
			return self.__diameter/2

	def getMass(self):
		return self.__density*self.getVolume()

	def getCOM(self):
		if self.__noseType == noseTypes[0]:
			V = self.getVolume()
			d = self.__thickness
			R2 = self.__diameter/2  # Outer radius of cone
			H2 = self.__height
			H1 = H2 - d
			R1 = R2*H1/H2
			a = 1/2*H2**2 - 2*H2**3/(3*H1) + H2**4/(4*H1**2)
			return np.pi/V*((H2*R2)**2/12 - a*R1**2)  # COM relative to bottom surface of nose
		elif self.__noseType == noseTypes[1]:
			d = self.__thickness
			R2 = self.__diameter/2
			R1 = R2 - d
			return 3/8*(R2**4 - R1**4)/(R2**3 - R1**3)  # COM relative to bottom surface of nose

	@staticmethod
	def from_file(file, noseType):
		if noseType.lower() == noseTypes[0]:  # Conic
			diameter = find_parameter(file, "diameter")
			height = find_parameter(file, "height")
			density = find_parameter(file, "density")
			thickness = find_parameter(file, "thickness")
			return Nose(eval(noseType.lower()), eval(diameter), eval(height), eval(thickness), eval(density))
		elif noseType.lower() == noseTypes[1]:
			diameter = find_parameter(file, "diameter")
			thickness = find_parameter(file, "thickness")
			density = find_parameter(file, "density")
			return Nose(eval(noseType.lower()), eval(diameter), eval(thickness), eval(density))
		else:
			print("ERROR: invalid nose type encountered at initialization. Please check your spelling.")
			exit(1)


class Body:

	def __init__(self, diameter, length, thickness, density):
		self.__diameter = diameter
		self.__length = length
		self.__density = density
		self.__thickness = thickness
		print("Body initialized!")

	def __str__(self):
		D = str(self.__diameter)
		l = str(self.__length)
		d = str(self.__thickness*1e3)
		m = str(self.getMass())
		rho = str(self.__density)
		return "Diameter: " + D + " m\n" + "Length: " + l + " m\n" + "Thickness: " + d + " mm\n" + "Mass: " + m + \
			   " kg\n" + "Density: " + rho + " kgm^-3 "

	# Member functions
	def getVolume(self):
		d = self.__thickness
		R2 = self.__diameter/2
		R1 = R2 - d
		h = self.__length
		return np.pi*(R2**2 - R1**2)*h

	def getCavityVolum(self):
		d = self.__thickness
		R = self.__diameter/2 - d
		l = self.__length
		return np.pi*l*R**2

	def getLength(self):
		return self.__length

	def getDiameter(self):
		return self.__diameter

	def getMass(self):
		return self.__density*self.getVolume()

	def getCOM(self):
		l = self.__length
		return -l/2  # COM relative to top of body

	@staticmethod
	def from_file(file):
		diameter = find_parameter(file, "diameter")
		length = find_parameter(file, "length")
		density = find_parameter(file, "density")
		thickness = find_parameter(file, "thickness")
		return Body(eval(diameter), eval(length), eval(thickness), eval(density))


class Fin:

	def __init__(self, *args):
		self.__semiChord = args[0]  # Fin semi span
		self.__rootChord = args[1]  # Fin root chord
		self.__tipChord = args[2]  # Fin tip chord
		self.__angle = args[3]  # Angle of ray from body to top outer edge
		self.__thickness = args[4]
		self.__density = args[5]
		print("Fin initialized!")

	def __str__(self):
		Chord = str(self.__semiChord)
		d = str(self.__thickness*1e3)
		m = str(self.getMass())
		rho = str(self.__density)
		return "Semi chord: " + Chord + " m\n" + "Thickness: " + d + " mm\n" + "Mass: " + m + " kg\n" \
			   + "Density: " + rho + " kgm^-3"

	# Member functions

	def getVolume(self):
		l1 = self.__rootChord
		l2 = self.__tipChord
		cord = self.__semiChord
		d = self.__thickness
		return 1/2*(l1 + l2)*cord*d

	def getLength(self):
		return self.__rootChord

	def getMass(self):
		return self.__density*self.getVolume()

	def getCOM(self):
		cord = self.__semiChord
		l1 = self.__rootChord
		l2 = self.__tipChord
		a = self.__angle
		return -(cord/np.tan(a) + 3/4*(l2 - l1))  # COM relative to top edge of fin

	def getCOP(self):
		lr = self.__rootChord
		lt = self.__tipChord

	@staticmethod
	def from_file(file):
		cord = find_parameter(file, "cord")
		length1 = find_parameter(file, "length1")
		length2 = find_parameter(file, "length2")
		angle = find_parameter(file, "angle")
		density = find_parameter(file, "density")
		thickness = find_parameter(file, "thickness")
		return Fin(eval(cord), eval(length1), eval(length2), eval(angle), eval(thickness), eval(density))


class Motor:

	def __init__(self, *args):
		self.__name = args[0]
		self.__thrustMatrix = args[1]
		self.__timeArray = self.__thrustMatrix[:, 0]  # Assuming time values along 1st column
		self.__thrustArray = self.__thrustMatrix[:, 1]  # Assuming thrust values along 2nd column
		self.__diameter = args[2]
		self.__length = args[3]
		self.__initialPropellantMass = args[4]
		self.__frameMass = args[5]
		print("Interpolating thrust data...")
		self.__thrustFunction = interp1d(self.__timeArray, self.__thrustArray)  # Linear Interpolation for thrust curve
		print("Calculating total impulse...")
		self.__totalImpulse = round(simps(self.__thrustArray, self.__timeArray), 3)
		print("\tTotal impulse: %1.2f Ns" % self.__totalImpulse)
		self.__exhaustSpeed = self.__totalImpulse/self.__initialPropellantMass
		self.__burnTime = self.__timeArray[-1]
		self.__avgThrust = round(self.__totalImpulse/self.__burnTime, 3)
		dt = self.__burnTime/1e4  # This is usually in order of milli seconds
		timeList = np.arange(self.__timeArray[0], self.__burnTime, dt)
		iterations = len(timeList)
		propellantMass = self.__initialPropellantMass
		massFlow = self.__thrustFunction(timeList)/self.__exhaustSpeed
		self.__propellantMassList = np.zeros(iterations)
		self.__propellantMassList[0] = propellantMass
		print("Calculating mass loss over the burn time of %1.2f s..."%self.__burnTime)
		for i in range(iterations - 1):
			propellantMass -= dt/2*(massFlow[i] + massFlow[i + 1])  # Trapezoid rule for integration of mdot over time
			self.__propellantMassList[i + 1] = propellantMass
		self.__propellantMassFunction = interp1d(timeList, self.__propellantMassList)
		print("Motor %s initialized!\n" % self.__name)

	def __str__(self):
		I = str(self.__totalImpulse/1e3)
		avg = str(self.__avgThrust)
		timeMax, Tmax = self.getMaxThrust()
		bTime = self.__burnTime
		name = self.getName()
		sep = (len(name) + 7)*'_'
		return "Motor: " + name + "\n" + sep + "\nTotal impulse: " + I + " kNs\n" + "Average thrust: " + avg + " kN" + \
			   "Maximum thrust: " + str(Tmax/1e3) + " kN,\tat time " + timeMax + " s\n" + "Burntime: " + bTime + "s \n"

	# Get functions
	def getName(self):
		return self.__name

	def getMass(self, t):
		return self.__propellantMassFunction(t) + self.__frameMass

	def getCOM(self, t):
		propM = self.__propellantMassFunction(t)
		frameM = self.__frameMass
		Mtot = frameM + propM
		l = self.__length
		return -1/Mtot*(frameM + propM**2/Mtot)*l/2  # COM relative to top of motor

	def getLength(self):
		return self.__length

	def getMaxThrust(self):
		i = np.argmax(self.__thrustArray)
		maxThrust = self.__thrustMatrix[i]
		return maxThrust[0], maxThrust[1]  # 1 timeAtMax, 2 maxThrust

	def getAverageThrust(self):
		return self.__avgThrust

	def getTotalImpulse(self):
		return self.__totalImpulse

	# Set functions

	# Auxiliary functions

	def thrust(self, t):
		if t <= self.__burnTime:
			return self.__thrustFunction(t)
		else:
			return 0

	def massFlow(self, t):
		return self.thrust(t)/self.__exhaustSpeed

	def plotPerformance(self):
		dt = self.__burnTime/1e4
		timeList = np.arange(self.__timeArray[0], self.__burnTime, dt)
		thrustArray = self.__thrustFunction(timeList)
		propellantMassArray = self.__propellantMassList
		# PLOT FORCE
		plt.suptitle(self.getName())
		plt.subplot(211)
		plt.plot(timeList, thrustArray, label='Thrust', c='r', lw='2')
		plt.title('Fuel mass during burn phase')
		plt.ylabel('thrust [N]')
		plt.xlabel('time [s]')
		plt.grid()
		plt.legend(loc='best')
		# PLOT PROPELLANT MASS LOSS
		plt.subplot(212)
		plt.plot(timeList, propellantMassArray, label='propellant mass', c='b', lw='2')
		plt.title('Fuel mass during burn phase')
		plt.ylabel('mass [kg]')
		plt.xlabel('time [s]')
		plt.grid()
		plt.legend(loc='best')
		plt.subplots_adjust(hspace=0.7)

		plt.show()

	@staticmethod
	def from_motorFile(motorFile):
		"""
			Read a file with motor specs.
			ASSUMPTIONS: -
			:param motorFile: The file with motor specs in proper format
			:return: motor
			"""
		name = find_parameter(motorFile, "name")
		diameter = find_parameter(motorFile, "diameter")
		length = find_parameter(motorFile, "length")
		propMass = find_parameter(motorFile, "propellant_mass")
		frameMass = find_parameter(motorFile, "frame_mass")
		motorFile = open(motorFile, 'r')

		thrust = [[0, 0]]
		with motorFile as fp:
			for i, line in enumerate(fp):
				if i > 5:
					row = line.split()
					t = eval(row[0])
					f = eval(row[1])
					thrust.append([t, f])
		thrust = np.array(thrust)
		return Motor(name, thrust, eval(diameter), eval(length), eval(propMass), eval(frameMass))


class Payload:

	def __init__(self, mass, placement, name="Bertha"):
		self.__mass = mass
		# Assuming placement is COM relative to body top
		self.__placement = placement
		self.__name = name

	def getMass(self):
		return self.__mass

	def getCOM(self):
		return self.__placement  # Assuming placement is position of COM

	@staticmethod
	def from_file(file, name):
		mass = find_parameter(file, "mass")
		placement = find_parameter(file, "placement")
		return Payload(eval(mass), eval(placement), name)


# TODO: implement COP and Inertia tensor for each part
# TODO: Placement of relevant rocket parts can be solved by specifying placement in rocket file


class RocketSimple:

	def __init__(self, nose, body, fin, motor, payload, partsPlacement):
		self.__rocketParts = np.array([nose, body, fin, payload])
		self.__rocketMotor = motor
		self.__massOfRocketParts = np.array([part.getMass() for part in self.__rocketParts])
		self.__rocketMass = self.__massOfRocketParts.sum()
		self.__motorMass = self.__rocketMotor.getMass(0)
		# TOTAL MASS
		self.__mass = self.__rocketMass + self.__motorMass
		self.__noseCOM = nose.getCOM() - nose.getLength()
		self.__bodyCOM = body.getCOM() - body.getLength() - nose.getLength()
		# Assuming placement of fin is position of top edge relative to body top
		self.__finCOM = partsPlacement[0] + fin.getCOM() - nose.getLength()
		# Assuming placement of motor is position of motor top relative to body top
		self.__motorCOM = partsPlacement[1] + motor.getCOM(0) - nose.getLength()
		self.__payloadCOM = partsPlacement[2]
		self.__COMofRocketParts = np.array([self.__noseCOM, self.__bodyCOM, self.__finCOM,
											self.__payloadCOM])
		self.__rocketCOM = (self.__massOfRocketParts*self.__COMofRocketParts).sum()/self.__rocketMass
		# FINAL COM OF ROCKET
		self.__COM = (self.__rocketCOM*self.__rocketMass + self.__motorCOM*self.__motorMass)/self.__mass

	def getNose(self):
		return self.__rocketParts[0]

	def getBody(self):
		return self.__rocketParts[1]

	def getFin(self):
		return self.__rocketParts[2]

	def getMotor(self):
		return self.__rocketMotor

	def getPayload(self):
		return self.__rocketParts[3]

	def getMass(self, t):
		self.__mass = self.__rocketMass + self.__rocketMotor.getMass(t)
		return self.__mass

	def getCOM(self, t):
		mass = self.getMass(t)
		motorCOM = self.getMotor().getCOM(t)
		motorMass = self.getMotor().getMass(t)
		self.__COM = (self.__rocketCOM*self.__rocketMass + motorCOM*motorMass)/mass
		return self.__COM

	@staticmethod
	def from_file(file):
		pass


def find_parameter(file, parameter):
	File = open(file, 'r')
	arr = ["", ""]
	while arr[0] != parameter.lower():
		base = File.readline()
		if base == '':
			print("ERROR: Could not find parameter '" + parameter + "' in '" + file + "'.")
			exit(1)
		base = base.replace(" ", "")
		base = base.replace("\n", "")
		arr = base.split("=")
	File.close()
	return arr[1]

def read_report(file, c):
	print("Reading force report...")
	with open(file) as file:
		filestring = ""
		for i in range(3):
			next(file)
		for line in file:
			filestring += line
		filestring = filestring.replace('\n', ' ')
		filelist = filestring.split(" ")
		filelist.remove('')
		for i in range(len(filelist)):
			filelist[i] = float(filelist[i])
		reportdata = np.array(filelist).reshape((-1, 3)).transpose()
		reportdata[0] *= c
		return reportdata

def plot_force_report(report):
	print("plotting force report...")
	fig, (axs1, axs2) = plt.subplots(2, sharex = True)
	fig.suptitle("Aerodynamic forces / AoA")
	axs1.plot(report[0] , report[1], color = "red", label = "Lift")
	axs2.plot(report[0], report[2], color = "blue", label = "Drag")
	axs1.legend(loc = "upper left")
	axs2.legend(loc = "lower left")
	axs1.grid()
	axs2.grid()
	axs1.set_title("Lift force")
	axs2.set_title("Drag force")
	axs1.set_ylabel("Force [N]")
	axs2.set_ylabel("Force [N]")
	plt.show()

motor1 = Motor.from_motorFile('motor.dot')
motor1.plotPerformance()
motor2 = Motor.from_motorFile('CesaroniM1520.dot')
motor2.plotPerformance()
motor3 = Motor.from_motorFile('CesaroniM1800.dot')
motor3.plotPerformance()
plot_force_report(read_report("sample_force_report.out", 0.5))
