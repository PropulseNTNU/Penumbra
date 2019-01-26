"""
Module containing models of the forces acting on a rocket

Version: 1.0
Last edit: 17.11.2018

--Propulse NTNU--
"""
import numpy as np
from scipy.constants import R, g, atmosphere

T0 = 20 + 273 # Temperature at sea level [K]
P0 = atmosphere # Air pressure at sea level [Pa]
m = 29e-3  # Molecular mass of Air [kg]
rho0 = P0/(R*T0/m)  # Air density at sea level [kg*m^-3]
h = R*T0/(m*g)  # Height constant of Air ~ 1e4 [m]

# Forces

def SAMdrag(rocket, position, linearVelocity):
    """
        Assumptions:- AoA ~ 0
                    - Quadratic drag F ~ -kv^2
    :param rocket: [rocket class] The rocket object
    :param position: [np.array] The position vector in world coordinates
    :param linearVelocity: [np.array] The current rocket velocity in world coord. (with wind)
    :return: [np.array] drag force in the world frame
    """
    z = abs(position[2])  # Vertical position of rocket
    Cd = rocket.getCd()
    Aref = np.pi*(rocket.getBody().getDiameter()/2)**2
    k = 1/2*rho0*Aref*Cd*np.exp(-z/h)

    return -k*np.linalg.norm(linearVelocity)*linearVelocity

def SAMlift(rocket, position, linearVelocity, AoA):
    """
    :param rocket: [rocket class] The rocket object
    :param position: [np.array] The position vector in world coordinates
    :param linearVelocity: [np.array] The current rocket velocity in world coord. (with wind)
    :return: [np.array] Lift force in the body frame
    """
    z = abs(position[2])  # Vertical position of rocket
    Cn = rocket.getCn(AoA)
    Aref = np.pi*(rocket.getBody().getDiameter()/2)**2
    k = 1/2*rho0*Aref*Cn*np.exp(-z/h)
    speed = np.linalg.norm(linearVelocity)

    return k*speed**2
