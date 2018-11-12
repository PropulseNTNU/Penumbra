"""
This is a module that contains functions for file reading.

Last edit: 11.11.2018
"""
from numpy import linspace, reshape, flip, loadtxt, sin, pi


def find_parameter(file, parameter):
	File = open(file, 'r')
	arr = ["", ""]
	while arr[0] != parameter.lower():
		base = File.readline()
		if base == '':
			print("ERROR: Could not find parameter '" + parameter + "' in '" + file + "'.")
			return False
		base = base.replace(" ", "")
		base = base.replace("\n", "")
		arr = base.split("=")
	File.close()
	return arr[1]


def unwrap_report(file, T, alpha_max, delta_v, v0):
	f = 1/T
	report = loadtxt(file, skiprows = 7).transpose()[1:]
	c = len(report.transpose())
	c -= c % T

	def make2d(report, c, T, idx):
		arr = report[idx]
		arr = arr[:c]
		arr = reshape(arr, (-1, T))
		for i in range(len(arr)):
			if i % 2 == 0:
				arr[i] = flip(arr[i])
		return arr.transpose()

	drag = make2d(report, c, T, 0)
	lift = make2d(report, c, T, 1)
	mom_x = make2d(report, c, T, 2)
	mom_y = make2d(report, c, T, 3)
	mom_z = make2d(report, c, T, 4)

	# Speed axis seems to be working fine regardless of T
	# Alpha axis acts weird when given other integers than 5
	alpha_axis = alpha_max * sin(linspace(0, T, T) * pi * f/2)
	speed_axis =  v0 + (linspace(0, c / T - 1, c // T) * delta_v)

	return drag, lift, mom_x, mom_y, mom_z, alpha_axis, speed_axis
