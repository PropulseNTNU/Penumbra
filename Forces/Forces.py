"""
Module containing models of the forces acting on a rocket

Version: 1.0
<<<<<<< b3d8334e28cbde7ced206b8c6d6bfc6785e68406
Last edit: 12.11.2018
=======
Last edit: 01.11.2018
>>>>>>> Restructured the folders, the implementation of the CFD rocket has started.

--Propulse NTNU--
"""
import numpy as np
from scipy.constants import R, g

<<<<<<< b3d8334e28cbde7ced206b8c6d6bfc6785e68406
#TODO Implement method to initiate T0 and P0
# Constants
#T0 = find_parameter("environment.dot", "temperature") + 273  # Temperature of Air [K]
#P0 = find_parameter("environment.dot", "pressure")  # Barometric pressure at sea level [Pa]
=======
# Constants
T0 = find_parameter("environment.dot", "temperature") + 273  # Temperature of Air [K]
P0 = find_parameter("environment.dot", "pressure")  # Barometric pressure at sea level [Pa]
>>>>>>> Restructured the folders, the implementation of the CFD rocket has started.
m = 29e-3  # Molecular mass of Air [kg]
#rho0 = P0/(R*T0/m)  # Air density at sea level [kg*m^-3]
#h = R*T0/(m*g)  # Height constant of Air ~ 1e4 [m]


# Hei detter er ny esidfj if



# Forces

def gravity(rocket, t):
	"""
	Gravitational force on the rocket
	Assumptions:- Homogeneous gravity field in -z direction
				- g to be a constant (very good approx. up to 3000 m, off by 0.1%)

	:param rocket: [rocket class] The rocket object
	:param t: [float] point in time [s]
	:return: [np.array] gravity vector in the inertial frame [N].
	"""
	return np.array([0, 0, -rocket.getMass(t)*g])


def thrust(rocket, t):
	"""
	:param rocket: [rocket class] The rocket object
	:param t: [float] point in time [s]
	:return: [np.array] thrust force at time t in rocket frame
	"""
	Thrust = rocket.getMotor().thrust(t)
	return np.array([Thrust, 0, 0])


def SAMdrag(rocket, position, linearVelocity):
	"""
		Assumptions:- AoA ~ 0
					- Quadratic drag F ~ -kv^2
	:param rocket: [rocket class] The rocket object
	:param position: [np.array] The position vector in world coordinates
	:param linearVelocity: [np.array] The current linear velocity
	:return: [np.array] The drag force (SAM) in the world frame
	"""
	z = position[2]  # Vertical position of rocket
	Cd = rocket.getCd()
	Aref = np.pi*(rocket.getBody().getDiameter()/2)**2
	k = 1/2*rho0*Aref*Cd*np.exp(-z/h)

	return -k*np.linalg.norm(linearVelocity)*linearVelocity