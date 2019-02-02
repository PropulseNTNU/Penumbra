import sys
sys.path.append('../Rocket/')
sys.path.append('../Forces/')
sys.path.append('../Trajectory/')
import numpy as np
import math
import scipy.linalg as splinalg
import scipy.integrate as spintegrate
import Kinematics
import Forces

epsilon = 1e-10

def calculateTrajectory(rocket, initialInclination, launchRampLength, timeStep, simulationTime):
    # x is the state of the vector
    # x = [position, quaternion, linear velocity, angular velocity]
    (x0, initialDirection) = initialState(rocket, initialInclination)
    t, x = integrateEquationsMotion(rocket, x0, launchRampLength, initialDirection, timeStep, simulationTime)
    (position, euler, linearVelocity, angularVelocity) = unwrapState(x)
     # Calculate AoA during flight, transform velocity to world frame for plot
    AoA = np.zeros(len(t))
    velocity = np.zeros(shape=(len(t),3)) # in world frame
    for i in range(len(t)):
        RotationBody2Inertial = Kinematics.Ryzx(euler[i][0], euler[i][1], euler[i][2])
        xAxisBody = RotationBody2Inertial[:,0]
        velocity[i] = RotationBody2Inertial @ linearVelocity[i]
        dirWindVelocity = (velocity[i]/(np.linalg.norm(velocity[i]) + epsilon))
        aoa = np.arccos(np.dot(dirWindVelocity, xAxisBody))
        AoA[i] = aoa
    # thrust, gravity, drag, lift
    return t, position, euler, AoA, velocity, angularVelocity

def initialState(rocket, initialInclination):
    initialPitch = np.pi/2-initialInclination
    R = Kinematics.Ryzx(initialPitch, 0, 0)
    initialDirection= R[:,0]
    initialPosition = rocket.getLength()*initialDirection
    initialQuaternion = Kinematics.euler2quaternion(initialPitch, 0, 0)
    initialLinearVelocity = np.array([0, 0, 0])
    initialAngularVelocity = np.array([0, 0, 0])
    x0 = np.concatenate((initialPosition, initialQuaternion,
    initialLinearVelocity, initialAngularVelocity))
    return (x0, initialDirection)

def integrateEquationsMotion(rocket, x0, launchRampLength, initialDirection, timeStep, simulationTime):
    N = math.ceil(simulationTime/timeStep)
    t = np.arange(0, simulationTime + timeStep, timeStep)
    x = np.zeros(shape=(N,len(x0)))
    sol = RK4(equationsMotion, 0, simulationTime, timeStep, x0, RHS_args=(rocket, launchRampLength, initialDirection))
    return t, sol

def equationsMotion(x, t, rocket, launchRampLength, initialDirection):
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
    dPosition = RotationBody2Inertial @ linearVelocity.T
    dQuaternion = Kinematics.quaternionGradient(quaternion) @ angularVelocity.T
    # forces in the body frame
    thrust = np.array([rocket.getMotor().thrust(t), 0, 0])
    gravity = RotationInertial2Body @ np.array([0, 0, rocket.getMass(t)*Forces.g])
    # aerodynamic forces
    windVelocity = np.array([0, 0, 0])
    # Add wind to current rocket velocity to get total air velocity
    airVelocity = dPosition + windVelocity
    airSpeed = np.linalg.norm(airVelocity)
    xAxisBody = RotationBody2Inertial[:,0]
    dirWindVelocity = (airVelocity/(np.linalg.norm(airVelocity) + epsilon))
    AoA = np.arccos(np.dot(dirWindVelocity, xAxisBody))
    dirDragBody = RotationInertial2Body @ (-dirWindVelocity.T)
    projectedDragBody = np.array([0, dirDragBody[1], dirDragBody[2]])
    dirProjectedDragBody = projectedDragBody/(np.linalg.norm(projectedDragBody) + epsilon)
    dirLiftBody = np.sin(AoA)*np.array([1, 0, 0]) + np.cos(AoA)*dirProjectedDragBody
    aeroForces = rocket.getAeroForces(AoA, position, airVelocity)
    drag = RotationInertial2Body @ aeroForces[0].T
    lift = aeroForces[1]*dirLiftBody*0
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
    totalForce = thrust + gravity + drag + lift
    if stillAtLaunchRamp:
        totalForce = np.array([totalForce[0], 0, 0])
        totalMoment = np.array([0, 0, 0])
    else:
        arm = rocket.getCOP(AoA) - rocket.getCOM(t)
        totalMoment = np.cross(arm, drag + lift)
    genForceBody = H.T @ np.concatenate((totalForce, totalMoment))
    # find dx
    genVelocity = np.concatenate((linearVelocity, angularVelocity))
    rhs = genForceBody - CBody @ genVelocity.T
    dGeneralizedVelocity = np.linalg.solve(IBody, rhs)
    dx = np.concatenate((dPosition, dQuaternion, dGeneralizedVelocity))
    # thrust, gravity, drag, lift
    return dx

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

    # Find times to evaluate and initialize array of positions
    timelist = np.arange(tmin, tmax + dt, dt)
    stateMatrix = np.zeros((len(timelist), len(w0)))

    # Initialize
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
        s1 = g(w, t)
        s2 = g(w + dt / 2 * s1, t + dt / 2)
        s3 = g(w + dt / 2 * s2, t + dt / 2)
        s4 = g(w + dt * s3, t + dt)

        w = w + dt / 6 * (s1 + 2 * s2 + 2 * s3 + s4)
        stateMatrix[i] = w
        #print("Iteration %5d/%d" % (i, steps - 1))

    # Return state at every instance (one row is one instance)
    return stateMatrix

def unwrapState(x):
    position = x[:, 0:3]
    quaternion = x[:, 3:7]
    linearVelocity = x[:, 7:10]
    angularVelocity = x[:, 10:13]
    euler = np.zeros(shape=(len(quaternion), 3))
    for index in range(len(quaternion)):
        euler[index,:] = Kinematics.quaternion2euler(quaternion[index,:])
    return (position, euler, linearVelocity, angularVelocity)
