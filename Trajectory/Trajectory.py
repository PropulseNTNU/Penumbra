"""
A library to solve the equations of motion for a 6 DOF rocket object.

2 frames of reference:
    - World frame, NED convention (North(x)-East(y)-Down(z)) is the inertial frame
    - Body frame, origin at nose-tip (x-axis is along the body (longitudinal),
      y and z completes the Right hand system)

State of rocket: [position, quaternion, linear Velocity, angular Velocity]

"""
import sys
sys.path.append('../Rocket/')
sys.path.append('../Forces/')
import numpy as np
import scipy.linalg as splinalg
import scipy.integrate as spintegrate
import Kinematics
import Forces

# Avoid division by 0 by adding epsilon to all denominators
epsilon = 1e-10

def calculateTrajectory(rocket, initialInclination, launchRampLength, timeStep, simulationTime):
    # x is the state of the rocket
    # x = [position, quaternion, linear velocity, angular velocity]
    (x0, initialDirection) = initialState(rocket, initialInclination)
    t, x, AoA, forces, aero_coeff = integrateEquationsMotion(rocket, x0, launchRampLength, initialDirection, timeStep, simulationTime)
    (position, euler, linearVelocity, angularVelocity) = unwrapState(x)
    n = len(t)
    drag = np.array([forces[i][:,0] for i in range(n)])
    lift = np.array([forces[i][:,1] for i in range(n)])
    gravity = np.array([forces[i][:,2] for i in range(n)])
    thrust = np.array([forces[i][:,3] for i in range(n)])
    # Transform velocity to world frame for plot
    velocity = np.zeros(shape=(len(t),3)) # in world frame
    for i in range(len(t)):
        RotationBody2Inertial = Kinematics.Ryzx(euler[i][0], euler[i][1], euler[i][2])
        velocity[i] = RotationBody2Inertial @ linearVelocity[i]

    return t, position, euler, AoA, velocity, angularVelocity, drag, lift, gravity, thrust, aero_coeff

def initialState(rocket, initialInclination):
    # Initial inclination of rocket
    initialPitch = np.pi/2-initialInclination
    # Calculate rotation matrix (the rocket orientation measured in world coords.)
    # at this initial inclination
    R = Kinematics.Ryzx(initialPitch, 0, 0)
    # rocket longitudinal axis in world frame (x-axis of body)
    initialDirection= R[:,0]
    initialPosition = rocket.getLength()*initialDirection
    initialQuaternion = Kinematics.euler2quaternion(initialPitch, 0, 0)
    initialLinearVelocity = np.array([0, 0, 0])
    initialAngularVelocity = np.array([0, 0, 0])
    # Initial state vector
    x0 = np.concatenate((initialPosition, initialQuaternion,
    initialLinearVelocity, initialAngularVelocity))
    return (x0, initialDirection)

def integrateEquationsMotion(rocket, x0, launchRampLength, initialDirection, timeStep, simulationTime):
    t = np.arange(0, simulationTime + timeStep, timeStep)
    x = np.zeros(shape=(len(t),len(x0)))
    sol, AoA, force, aero_coeff = RK4(equationsMotion, 0, simulationTime, timeStep, x0, RHS_args=(rocket, launchRampLength, initialDirection))
    return t, sol, AoA, force, aero_coeff

def equationsMotion(x, t, rocket, launchRampLength, initialDirection):
    """
    x: [np.array] the current state of rocket
    t: [float] at time t
    rocket: [rocket object] the rocket object
    launchRampLength: [float] The length of the launch ramp
    initialDirection: [np.array] the initial direction of rocket (in world coords)

    dx: Derivative of state x [np.array]
    AoA: Angle of attack in radians [float]
    forceMatrix: [3x4 matrix] A collection of the forces (one force vector for each column)

    return: dx, AoA, forceMatrix
    """
    position = x[0:3]
    quaternion = x[3:7]
    linearVelocity = x[7:10]
    angularVelocity = x[10:13]
    stillAtLaunchRamp = False
    # determine if whether at launch ramp or not
    if np.dot(position, initialDirection) <= launchRampLength + rocket.getLength():
        stillAtLaunchRamp = True
    else:
        stillAtLaunchRamp = False
    # dPosition and dQuaternion
    RotationBody2Inertial = Kinematics.Rquaternion(quaternion)
    RotationInertial2Body = RotationBody2Inertial.T
    # Velocity of rocket in world frame
    dPosition = RotationBody2Inertial @ linearVelocity.T
    dQuaternion = Kinematics.quaternionGradient(quaternion) @ angularVelocity.T
    # forces in the body frame
    thrust = np.array([rocket.getMotor().thrust(t), 0, 0])
    gravityWorld = np.array([0, 0, rocket.getMass(t)*Forces.g])
    gravityBody = RotationInertial2Body @ gravityWorld
    # aerodynamic forces
    # Add wind to current rocket velocity to get total air velocity
    airVelocity = dPosition
    airSpeed = np.linalg.norm(airVelocity)
    xAxisBody = RotationBody2Inertial[:,0]
    dirWindVelocity = (airVelocity/(np.linalg.norm(airVelocity) + epsilon))
    # definition of angle of attack
    AoA = np.arccos(np.dot(dirWindVelocity, xAxisBody))
    # unit vector that points in drag direction (body coords.)
    dirDragBody = RotationInertial2Body @ (-dirWindVelocity.T)
    projectedDragBody = np.array([0, dirDragBody[1], dirDragBody[2]])
    dirProjectedDragBody = projectedDragBody/(np.linalg.norm(projectedDragBody) + epsilon)
    # unit vector that points in lift direction (body coords.)
    dirLiftBody = np.sin(AoA)*np.array([1, 0, 0]) + np.cos(AoA)*dirProjectedDragBody
    # Get drag and lift for current state
    Forces.updateCd_2(rocket, position, linearVelocity, AoA, rocket.getCompressibilityState())
    aeroForces = rocket.getAeroForces(AoA, position, airVelocity)
    drag = RotationInertial2Body @ aeroForces[0].T
    lift = aeroForces[1]*dirLiftBody
    # inertia matrix and coriolis matrix for equations of motion
    # seen from origin of body frame, not from center of mass (See Fossen)
    H = Kinematics.TransformationMatrix(rocket.getCOM(t))
    m = rocket.getMass(t)
    I = rocket.getInertiaMatrix(t)
    IBody = H.T @ splinalg.block_diag(m,m,m,I) @ H
    S1 = Kinematics.CrossProductMatrix(m*angularVelocity)
    S2 = Kinematics.CrossProductMatrix(I @ angularVelocity.T)
    CBody = H.T @ splinalg.block_diag(S1, -S2) @ H
    # obtain generalized forces seen from origin of body frame
    totalForce = thrust + gravityBody + drag + lift
    forceMatrix = np.array([drag, lift, gravityBody, thrust]).T
    aero_coeff = np.array([rocket.getCd(), rocket.getCn(AoA)])
    # If rocket is on launch ramp, don't allow it to rotate, only accelerate (fixed to ramp)
    if stillAtLaunchRamp:
        totalForce = np.array([totalForce[0], 0, 0])
        totalMoment = np.array([0, 0, 0])
    else:
        # After launch ramp, allow it to rotate (now calculating torques about COM)
        arm = rocket.getCOP(AoA) - rocket.getCOM(t)
        totalMoment = np.cross(arm, drag + lift)
    genForceBody = H.T @ np.concatenate((totalForce, totalMoment))
    # find dx
    genVelocity = np.concatenate((linearVelocity, angularVelocity))
    rhs = genForceBody - CBody @ genVelocity.T
    dGeneralizedVelocity = np.linalg.solve(IBody, rhs)
    dx = np.concatenate((dPosition, dQuaternion, dGeneralizedVelocity))

    return dx, AoA, forceMatrix, aero_coeff

# Solving simultaneous diff. equations
def RK4(RHS, tmin, tmax, dt, w0, RHS_args=0):
    """
    Runge-Kutta ODE solver of order 4, solving the equation system given by RHS

    :param RHS: RHS(w, t); Derivative of the state w (right hand side of system)
    :param tmin: Float; starting time, >= 0
    :param tmax: Float; ending time, > tmin
    :param dt: Float; time step, > 0
    :param w0: RHS-parameter-type; the initial state of the system
    :param RHS_args: [tupple] if RHS has several arguments, insert in order as a tupple.
    :return: np.array[] ; matrix that contains the states at every instance (as row vectors)
    """

    # array of time, initialize arrays to store state vectors (one state for each row)
    timelist = np.arange(tmin, tmax + dt, dt)
    stateMatrix = np.zeros((len(timelist), len(w0)))
    forceMatrix = np.zeros((len(timelist), 3, 4))
    aero_coeffMatrix = np.zeros((len(timelist), 2))
    aoa = np.zeros(len(timelist))

    # Initialize with initial conditions.
    w = w0
    stateMatrix[0] = w
    steps = len(timelist)

    if RHS_args:
        def g(w, t): return RHS(w, t, *RHS_args)
    else:
        def g(w, t): return RHS(w, t)

    # Runge-Kutta algorithm
    for i in range(1, steps):
        t = timelist[i]
        s1, AoA, force, aero_coeff = g(w, t)
        s2 = g(w + dt / 2 * s1, t + dt / 2)[0] # get dx only
        s3 = g(w + dt / 2 * s2, t + dt / 2)[0] # get dx only
        s4 = g(w + dt * s3, t + dt)[0] # get dx only

        w = w + dt / 6 * (s1 + 2 * s2 + 2 * s3 + s4)
        stateMatrix[i] = w  # Store new state
        forceMatrix[i] = force  # Store forces
        aero_coeffMatrix[i] = aero_coeff # Store aero-coefficients
        aoa[i] = AoA  # store AoA

    # Return state, AoA, forces and aero_coefficients at every instance (np.arrays)
    return stateMatrix, aoa, forceMatrix, aero_coeffMatrix

def unwrapState(x):
    """
    A simple routine to unwrap the state vector to its costituent parts (position, orientation, velocity etc.)
    """
    position = x[:, 0:3]
    quaternion = x[:, 3:7]
    linearVelocity = x[:, 7:10]
    angularVelocity = x[:, 10:13]
    euler = np.zeros(shape=(len(quaternion), 3))
    for index in range(len(quaternion)):
        euler[index,:] = Kinematics.quaternion2euler(quaternion[index,:])
    return (position, euler, linearVelocity, angularVelocity)
