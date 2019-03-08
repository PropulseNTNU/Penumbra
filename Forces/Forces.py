"""
Module containing models of aero forces acting on a rocket

Version: WIP
Last edit: 04.02.2019

--Propulse NTNU--
"""
import numpy as np
from scipy.constants import R, g, atmosphere

T0 = 20 + 273 # Temperature at sea level [K]
P0 = atmosphere # Air pressure at sea level [Pa]
m = 29e-3  # Molecular mass of Air [kg]
rho0 = P0/(R*T0/m)  # Air density at sea level [kg/m^3]
h = R*T0/(m*g)  # Height constant of Air ~ 1e4 [m]
nu = 1.511e-5  # Kinematic viscosity of air [m^2/s]
c = 343  # Speed of sound (at 293K) [m/s]

# Forces
def Drag1(rocket, position, linearVelocityBody, AoA):
    """
    Reference: OpenRocket techDoc, section 3.4.2
       Assuming contribution to skin drag is component of velocity along body.
    :param rocket: [rocket class] The rocket object
    :param position: [np.array] The position vector in world coordinates
    :param linearVelocity: [np.array] The current rocket velocity in body coord. (with wind)
    :return: [np.array] drag force in the world frame [N]

    """
    z = abs(position[2])  # Vertical position of rocket
    velocity = np.array([linearVelocityBody[0], 0, 0])  # component along body x-axis that contributes
    speed = np.linalg.norm(velocity)
    M = speed/c
    AwetNose = rocket.getNose().getSurfaceArea()
    AwetBody = rocket.getBody().getSurfaceArea()
    N = rocket.getNumberOfFins()
    AwetFins = 2*N*rocket.getFin().getSurfaceArea()
    Awet = AwetNose + AwetBody + AwetFins
    D = rocket.getBody().getDiameter()
    R = speed*D/nu  # Reynold's number (Kinematic viscosity)
    Rcrit = 51*(100e-6/D)**(-1.039)
    Cf = 0
    # Conditions for different Reynold's number
    if R < 1e4:
        Cf = 1.48e-2
    elif 1e4 < R < Rcrit:
        Cf = 1/(1.5*np.log(R)-5.6)**2
    else:
        Cf = 0.032*(100e-6/D)**0.2

    # Conditions for different speeds (subsonic/supersonic)
    if M < 0.8:
        Cf = Cf*(1 - 0.1*M**2)
    else:
        Cf = Cf/(1 + 0.15*M**2)**0.58

    k = 1/2*rho0*Awet*Cf*np.exp(-z/h)

    return -k*speed*velocity

def Drag2(rocket, position, linearVelocityBody, AoA):
    """
    Reference: "Estimating the dynamic and aerodynamic paramters of
    passively controlled high power rockets for flight simulaton" by Simon B. .. Feb 2009

    :param rocket: [rocket class] The rocket object
    :param position: [np.array] The position vector in world coordinates
    :param linearVelocity: [np.array] The current rocket velocity in body coord. (with wind)
    :return: [float] Lift force in the body frame [N]
    """
    z = abs(position[2])  # Vertical position of rocket
    velocity = np.array([linearVelocityBody[0], 0, 0])  # component along body x-axis that contributes
    speed = np.linalg.norm(velocity)
    M = speed/c  # Mach number
    Rcrit = 5e5
    # For body and nose
    Lb = rocket.getBody().getLength() + rocket.getNose().getLength()
    R = speed*Lb/nu  # Reynold's number (Kinematic viscosity)
    B = Rcrit*(0.074/(R**0.2) - 1.328/(R**0.5))
    Cfb = 0
    # Conditions for different R
    if R <= Rcrit:
        Cfb = 1.328/(R**0.5)
    else:
        Cfb = 0.074/(R**0.2) - B/R

    # For fins
    Lf = rocket.getFin().getSemiChord()
    R = R*Lf/Lb
    B = Rcrit*(0.074/(R**0.2) - 1.328/(R**0.5))
    Cff = 0
    # Conditions for different R
    if R <= Rcrit:
        Cff = 1.328/(R**0.5)
    else:
        Cff = 0.074/(R**0.2) - B/R

    # Interference term (between body and fins)
    D = rocket.getBody().getDiameter()
    N = rocket.getNumberOfFins()
    Cdint = 2*Cff*(1 + 2)

    Cf = Cfb + Cff
    # TODO finish this
    return Cf

def SAMdrag(rocket, position, linearVelocityWorld):
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

    return -k*np.linalg.norm(linearVelocityWorld)*linearVelocityWorld

def SAMlift(rocket, position, linearVelocityWorld, AoA):
    """
    :param rocket: [rocket class] The rocket object
    :param position: [np.array] The position vector in world coordinates
    :param linearVelocity: [np.array] The current rocket velocity in world coord. (with wind)
    :return: [float] Lift force in the body frame [N]
    """
    z = abs(position[2])  # Vertical position of rocket
    Cn = rocket.getCn(AoA)  # Lift coefficient
    Aref = np.pi*(rocket.getBody().getDiameter()/2)**2
    k = 1/2*rho0*Aref*Cn*np.exp(-z/h)
    speed = np.linalg.norm(linearVelocityWorld)

    return k*speed**2
