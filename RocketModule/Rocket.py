"""
Rocket module - The definition of the rocket and its constituent parts

Version: WIP
Last edit: 17.10.2018

--Propulse NTNU--
"""
import numpy as np
import scipy.integrate as quadrature
import scipy.interpolate as interpolate


class NoseCone:

	def __init__(self, diameter, length, thickness, density):
		# Member variables
		self.__diameter = diameter
		self.__length = length
		self.__thickness = thickness
		self.__density = density

	def __str__(self):
		D = str(self.__diameter)
		l = str(self.__length)
		d = str(self.__thickness)
		m = str(self.getMass())
		rho = str(self.__density)
		return "Diameter: " + D + "m\n" + "Length: " + l + "m\n" + "Thickness: " + d + "m\n\n" + "Mass: " + m + "kg\n" \
			   + "Density: " + rho + "kgm^-3 "

	# Member functions
	def getVolume(self):
		R = self.__diameter/2
		h = self.__length
		d = self.__thickness
		return np.pi/3*(h*R**2 - (h - d)*(R - 2*d)**2)

	def getMass(self):
		return self.__density*self.getVolume()

	def getCOM(self):
		pass

	def getCOP(self):
		pass

	@staticmethod
	def from_file(file):
		diameter = find_parameter(file, "diameter")
		length = find_parameter(file, "length")
		density = find_parameter(file, "density")
		thickness = find_parameter(file, "thickness")
		return NoseCone(diameter, length, thickness, density)


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
		return "Diameter: " + D + "m\n" + "Length: " + l + "m\n" + "Thickness: " + d + "m\n\n" + "Mass: " + m + "kg\n" \
			   + "Density: " + rho + "kgm^-3 "

	# Member functions
	def getVolume(self):
		R = self.__diameter/2
		h = self.__length
		d = self.__thickness
		return np.pi*(R**2 - (R - 2*d)**2)*h

	def getMass(self):
		return self.__density*self.getVolume()

	def getCOM(self):
		pass

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

	def __init__(self, cord, length, thickness, density, type):
		# Assuming fins are triangular/circular etc.
		self.__cord = cord  # cord length is the span of the fin
		self.__length = length
		self.__density = density
		self.__thickness = thickness
		self.__type = type

	def __str__(self):
		C = str(self.__cord)
		d = str(self.__thickness)
		m = str(self.getMass())
		rho = str(self.__density)
		return "Radius: " + C + "m\n" + "Thickness: " + d + "m\n\n" + "Mass: " + m + "kg\n" \
			   + "Density: " + rho + "kgm^-3"

	# Member functions

	def getVolume(self):
		R = self.__cord/2
		d = self.__thickness
		return (np.pi*R**2)*d

	def getMass(self):
		return self.__density*self.getVolume()

	def getCOM(self):
		pass

	def getCOP(self):
		pass

	@staticmethod
	def from_file(file):
		cord = find_parameter(file, "cord")
		length = find_parameter(file, "length")
		density = find_parameter(file, "density")
		thickness = find_parameter(file, "thickness")
		type = find_parameter(file, "type")
		return Fin(cord, length, thickness, density, type)


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