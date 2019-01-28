import sys
sys.path.append('../Rocket/')
sys.path.append('../Forces/')
sys.path.append('../Trajectory/')
import numpy as np
import math
import scipy.linalg as splinalg
import scipy.integrate as spintegrate
import Kinematics
import Rocket2 as Rocket
from scipy.constants import g

epsilon = 1e-10

def calculateTrajectory(rocket, initialInclination, launchRampLength, timeStep, simulationTime):
    # x is the state of the vector
    # x = [position, quaternion, linear velocity, angular velocity]
    (x0, initialDirection) = initialState(rocket, initialInclination)
    (t, x, AoA, thrust, gravity, drag, lift) = integrateEquationsMotion(rocket, x0, launchRampLength, initialDirection, timeStep, simulationTime)
    (position, euler, linearVelocity, angularVelocity) = unwrapState(x)
    return ((t, position, euler, linearVelocity, angularVelocity, AoA, thrust, gravity, drag, lift))

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
    t = np.arange(0,N*timeStep,timeStep)
    x = np.zeros(shape=(N,len(x0)))
    AoA = np.zeros(shape=(N-1,1))
    thrust = np.zeros(shape=(N-1,3))
    gravity = np.zeros(shape=(N-1,3))
    drag = np.zeros(shape=(N-1,3))
    lift = np.zeros(shape=(N-1,3))
    for index in range(len(t)):
        if index==0:
            x_current = x0
            x[index,:] = x_current
        else:
            previousIndex = index-1
            t_current = t[previousIndex]
            (dx, AoA_current, thrust_current, gravity_current, drag_current, lift_current) \
            = equationsMotion(rocket, x_current, t_current, launchRampLength, initialDirection)
            x_current = x_current + timeStep*dx #Forward Euler!
            x[index,:] = x_current
            AoA[previousIndex] = AoA_current
            thrust[previousIndex] = thrust_current
            gravity[previousIndex] = gravity_current
            drag[previousIndex] = drag_current
            lift[previousIndex] = lift_current
    return (t, x, AoA, thrust, gravity, drag, lift)

def equationsMotion(rocket, x, t, launchRampLength, initialDirection):
    position = x[0:3]
    quaternion = x[3:7]
    linearVelocity = x[7:10]
    angularVelocity = x[10:13]
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
    gravity = RotationInertial2Body @ np.array([0, 0, rocket.getMass(t)*g])
    # aerodynamic forces
    windVelocity = np.array([0, 0, 0])
    # Add wind to current rocket velocity to get total air velocity
    airVelocity = linearVelocity + windVelocity
    airSpeed = np.linalg.norm(airVelocity)
    xAxisBody = RotationBody2Inertial[:,0]
    dirWindVelocity = (airVelocity/(np.linalg.norm(airVelocity) + epsilon))
    AoA = np.arccos(np.dot(dirWindVelocity, xAxisBody))
    dirDragBody = RotationInertial2Body @ dirWindVelocity
    #TODO Problem with drag and lift directions
    projectedDragBody = np.array([0, dirDragBody[1], dirDragBody[2]])
    dirProjectedDragBody = projectedDragBody/(np.linalg.norm(projectedDragBody) + epsilon)
    dirLiftBody = np.sin(AoA)*np.array([1, 0, 0]) + np.cos(AoA)*dirProjectedDragBody
    aeroForces = rocket.getAeroForces(AoA, position, airVelocity)
    drag = RotationInertial2Body @ aeroForces[0]
    lift = -aeroForces[1]*dirLiftBody*0
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
    # AoA, thrust, gravity, drag, lift
    return dx

def unwrapState(x):
    position = x[:,0:3]
    quaternion = x[:,3:7]
    linearVelocity = x[:,7:10]
    angularVelocity = x[:,10:13]
    euler = np.zeros(shape=(len(quaternion),3))
    for index in range(len(quaternion)):
        euler[index,:] = Kinematics.quaternion2euler(quaternion[index,:])
    return (position, euler, linearVelocity, angularVelocity)
