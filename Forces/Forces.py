"""
Module containing models of aero forces acting on a rocket

Version: WIP
Last edit: 10.03.2019

--Propulse NTNU--
"""
import numpy as np
from scipy.constants import R, g, atmosphere

T0 = 27 + 273 # Temperature at sea level [K]
P0 = atmosphere # Air pressure at sea level [Pa]
m = 29e-3  # Molar mass of Air [kg]
rho0 = P0/(R*T0/m)  # Air density at sea level [kg/m^3]
h = R*T0/(m*g)  # Height constant of Air ~ 1e4 [m]
nu = 1.511e-5  # Kinematic viscosity of air [m^2/s]
c = 343  # Speed of sound (at 293K) [m/s]

# Forces
def updateCd(rocket, position, linearVelocityBody, AoA, enable_compressibility=True, deviation=0):
    """
    Reference:
    -"Estimating the dynamic and aerodynamic paramters of
    passively controlled high power rockets for flight simulaton" by Simon B. .. Feb 2009
    -"OpenRocket techDoc, section 3.4.2"

    "The comparison shows a good agreement between experimental and modelled data
    for vaules 0deg < AoA < 15deg, K=1" (page 12 in document)

    This function evaluates the drag coefficient Cd at the given state and updates the rocket object.
    NOTE:
            - Using skin friction model from OpenRocket [section 3.4.2]

    :param rocket: [rocket class] The rocket object
    :param position: [np.array] The position vector in world coordinates
    :param linearVelocity: [np.array] The current rocket velocity in body coord. (with wind)
    :param AoA: [float] the current angle of attack
    :param enable_compressibility: [bool] Set True if Compressibility should be taken into account. Default value: False
    """
    speed = np.linalg.norm(linearVelocityBody)
    if speed < 0.01:
        rocket.setCd(0)
        return True
    # Temperature at altitude z (drops 6K per kilometer)
    T = T0 - 6*position[2]/1000
    # Mach number at current state (speed of sound goes as sqrt(temperature))
    M = speed/c*np.sqrt(T0/T)
    Rcrit = 5e5
    # For body and nose
    Ln = rocket.getNose().getLength()
    Lb = rocket.getBody().getLength()
    Ltot = Lb + Ln
    D = rocket.getBody().getDiameter()
    R = speed*Ltot/nu  # Reynold's number (Kinematic viscosity)
    Cfb = 0
    # Conditions for different Reynold's number
    if R < 1e4:
        Cfb = 1.48e-2
    elif 1e4 < R < Rcrit:
        Cfb = 1/(1.5*np.log(R)-5.6)**2
    else:
        Cfb = 0.032*(100e-6/D)**0.2
    # Body drag contrib. (assuming no boat tail) (eq 41)
    Cd_fb = (1 + 60/(Ltot/D)**3 + (2.5e-3)*Lb/D)*(2.7*Ln/D + 4*Lb/D)*Cfb*0.45
    # Base drag contrib. (due to boundary layer seperation, low pressure region)
    if Cd_fb != 0:
        Cd_b = 0.029/np.sqrt(Cd_fb)
    else:
        Cd_b = 0
    # For fins
    Lm = rocket.getFin().getMidChord()
    Lr = rocket.getFin().getRootChord()
    Tf = rocket.getFin().getThickness()
    N = rocket.getNumberOfFins()
    Afe = 2*N*rocket.getFin().getSurfaceArea()
    Abe = rocket.getBody().getSurfaceArea()
    Afp = Afe + 1/2*D*Lr
    R = R*Lm/Ltot
    Cff = 0
    # Conditions for different Reynold's number
    if R < 1e4:
        Cff = 1.48e-2
    elif 1e4 < R < Rcrit:
        Cff = 1/(1.5*np.log(R)-5.6)**2
    else:
        Cff = 0.032*(100e-6/D)**0.2
    # Drag contrib. (43)
    Cd_f = Cff*((1 + D/(2*Ltot))*Abe + (1 + 2*Tf/Lm)*Afe)/(np.pi*D**2)
    #Cd_f = 2*Cff*(1 + 2*Tf/Lm)*4*N*Afp/(np.pi*D**2)

    # Interference term (between body and fins) (44)
    Cd_int = 2*Cff*(1 + 2*Tf/Lm)*4*N*(1/2*D*Lr)/(np.pi*D**2)

    # AoA corrections (assuming AoA is below 10 deg)
    # From body
    delta = 0.9 # experimental values (see page 14 in document referenced above)
    eta = 0.7
    # Fin section ratio Rs (total fin span / Body diameter)
    Rs = 2*np.sin(2*np.pi/N)*rocket.getTransversalExtension()/D
    kfb = 0.8065*Rs**2 + 1.1553*Rs
    kbf = 0.1935*Rs**2 + 0.8174*Rs + 1
    Cd_bA = 2*delta*AoA**2 + 3.6*eta*(1.36*Ltot - 0.55*Ln)/(np.pi*D)*(AoA**3)
    Cd_fA = (AoA**2)*4/(np.pi*D**2)*(1.2*Afp + 3.12*(kfb + kbf - 1)*Afe)

    # total zero AoA drag coefficient (eq 48)
    CD_0 = Cd_fb + Cd_b + Cd_f + Cd_int
    # FINAL CD
    CD = CD_0 + Cd_bA + Cd_fA + deviation

    # Compressibility:
    if enable_compressibility:
        if M < 0.90:
            CD /= np.sqrt(1 - M**2)
        elif 0.90 <= M < 1.1:
            CD /= np.sqrt(1 - (0.90)**2)
        else:
            CD /= np.sqrt(M**2 - 1)

    # Update Cd of rcoket object
    rocket.setCd(CD)

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
