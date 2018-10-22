"""
Rocket module - The definition of the rocket and its constituent parts

Version: WIP
Last edit: 21.10.2018

--Propulse NTNU--
"""
import numpy as np
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

	def __str__(self):

		if self.__noseType == noseTypes[0]:  # Conic
			D = str(self.__diameter)
			h = str(self.__height)
			d = str(self.__thickness*1000)
			m = str(round(self.getMass(), 2))
			rho = str(self.__density)
			return "Diameter: " + D + " m\n" + "Height: " + h + " m\n" + "Thickness: " + d + " mm\n" + "Mass: " + m + \
				   " kg\n" + "Density: " + rho + " kgm^-3"

		elif self.__noseType == noseTypes[1]:  # Dome
			D = str(self.__diameter)
			d = str(self.__thickness*1000)
			m = str(round(self.getMass(), 2))
			rho = str(self.__density)
			return "Diameter: " + D + " m\n" + " m\n" + "Thickness: " + d + " mm\n\n" + "Mass: " + m + \
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
			return np.pi/V*((H2*R2)**2/12 - a*R1**2)
		elif self.__noseType == noseTypes[1]:
			d = self.__thickness
			R2 = self.__diameter/2
			R1 = R2 - d
			return 3/8*(R2**4 - R1**4)/(R2**3 - R1**3)

	def getCOP(self):
		pass

	@staticmethod
	def from_file(file, noseType):
		if noseType.lower() == noseTypes[0]:  # Conic
			diameter = find_parameter(file, "diameter")
			height = find_parameter(file, "height")
			density = find_parameter(file, "density")
			thickness = find_parameter(file, "thickness")
			return Nose(noseType.lower(), diameter, height, thickness, density)
		elif noseType.lower() == noseTypes[1]:
			diameter = find_parameter(file, "diameter")
			thickness = find_parameter(file, "thickness")
			density = find_parameter(file, "density")
			return Nose(noseType.lower(), diameter, thickness, density)
		else:
			print("ERROR: invalid nose type encountered at initialization. Please check your spelling.")
			exit(1)


class Body:

	def __init__(self, diameter, length, thickness, density):
		self.__diameter = diameter
		self.__length = length
		self.__density = density
		self.__thickness = thickness

	def __str__(self):
		D = str(self.__diameter)
		l = str(self.__length)
		d = str(self.__thickness)
		m = str(self.getMass())
		rho = str(self.__density)
		return "Diameter: " + D + " m\n" + "Length: " + l + " m\n" + "Thickness: " + d + " m\n" + "Mass: " + m +\
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

	def getMass(self):
		return self.__density*self.getVolume()

	def getCOM(self):
		l = self.__length
		return l/2

	def getCOP(self):
		pass

	@staticmethod
	def from_file(file):
		diameter = find_parameter(file, "diameter")
		length = find_parameter(file, "length")
		density = find_parameter(file, "density")
		thickness = find_parameter(file, "thickness")
		return Body(diameter, length, thickness, density)


class Fin:

	def __init__(self, *args):
		self.__cord = args[0]
		self.__length1 = args[1]  # Length of fin along body
		self.__length2 = args[2]  # Length of fin along outer edge
		self.__angle = args[3]  # Angle of ray from body to top outer edge
		self.__thickness = args[4]
		self.__density = args[5]

	def __str__(self):
		Cord = str(self.__cord)
		d = str(self.__thickness)
		m = str(self.getMass())
		rho = str(self.__density)
		return "Cord: " + Cord + " m\n" + "Thickness: " + d + " m\n" + "Mass: " + m + " kg\n" \
			   + "Density: " + rho + " kgm^-3"

	# Member functions

	def getVolume(self):
		l1 = self.__length1
		l2 = self.__length2
		cord = self.__cord
		d = self.__thickness
		return 1/2*(l1 + l2)*cord*d

	def getMass(self):
		return self.__density*self.getVolume()

	def getCOM(self):
		cord = self.__cord
		l1 = self.__length1
		l2 = self.__length2
		a = self.__angle
		return cord/np.tan(a) + 3/4*(l2-l1)

	def getCOP(self):
		pass

	@staticmethod
	def from_file(file):
		cord = find_parameter(file, "cord")
		length1 = find_parameter(file, "length1")
		length2 = find_parameter(file, "length2")
		angle = find_parameter(file, "angle")
		density = find_parameter(file, "density")
		thickness = find_parameter(file, "thickness")
		return Fin(cord, length1, length2, angle, thickness, density)


class Motor:

	def __init__(self, structMass, fuelMass, thrustData):
		pass


# TODO: implement fins (triangular, circular, square), motor, recovery
# TODO: implement COM for each part
# TODO: implement COP for each part


class RocketSimple:

	def __init__(self, body, fin, nosecone, motor, recovery, payload, slug):
		self.__body = body
		self.__fin = fin
		self.__nosecone = nosecone
		self.__motor = motor
		self.__recovery = recovery
		self.__payload = payload
		self.__slug = slug


def find_parameter(file, parameter):
	File = open(file, 'r')
	arr = ["", ""]
	while arr[0] != parameter.lower():
		base = File.readline()
		if base == '':
			print("ERROR: Could not find parameter '" + parameter + "' in '" + file + "'.")
			exit(1)
		base = base.replace(" ", "")
		arr = base.split("=")
	File.close()
	return eval(arr[1])

nose = Nose.from_file('test.dot', 'conic')
print(nose.getCOM())
print(nose)