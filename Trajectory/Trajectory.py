import sys
sys.path.append('../Rocket/')
sys.path.append('../Forces/')
sys.path.append('../Trajectory/')
import numpy as np
import scipy.linalg as splinalg
import scipy.integrate as spintegrate
import Kinematics
import Rocket2 as Rocket
import Forces as Forces

def calculateTrajectory(rocket, initialInclination, launchRampLength, timeStep, simulationTime):
    # x is the state of the vector
    # x = [position, quaternion, linear velocity, angular velocity]
    (x0, initialDirection) = initialState(rocket, initialInclination)
    (t, x, AoA, thrust, gravity, drag, lift) = integrateEquationsMotion(rocket, x0, launchRampLength, initialDirection, timeStep, simulationTime)
    (position, quaternion, linearVelocity, angularVelocity) = unwrapState(x)
    return ((t, position, quaternion, linearVelocity, angularVelocity, AoA, thrust, gravity, drag, lift))

def initialState(rocket, initialInclination):
    R = Kinematics.Ryzx(initialInclination-90, 0, 0)
    initialDirection= R[:,0]
    initialPosition = rocket.getLength()*initialDirection
    initialQuaternion = Kinematics.euler2quaternion(initialInclination,0,0)
    initialLinearVelocity = np.array([0, 0, 0])
    initialAngularVelocity = np.array([0, 0, 0])
    x0 = np.concatenate((initialPosition, initialQuaternion,
    initialLinearVelocity, initialAngularVelocity))
    return (x0, initialDirection)

def integrateEquationsMotion(rocket, x0, launchRampLength, initialDirection, timeStep, simulationTime):
    N = ceil(simulationTime/timeStep)
    t = np.arange(0,N*timeStep,timeStep)
    x = np.zeros(shape=(N,len(x_0)))
    AoA = np.zeros(shape=(N-1,1))
    thrust = np.zeros(shape=(N-1,3))
    gravity = np.zeros(shape=(N-1,3))
    drag = np.zeros(shape=(N-1,3))
    lift = np.zeros(shape=(N-1,3))
    for index in range(len(time)):
        if index==0:
            x_current = x0
            x[index,:] = x_current
        else:
            previousIndex = index-1
            t_current = t[previousIndex]
            (dx, AoA_current, thrust_current, gravity_current, drag_current, lift_current) = equationsMotion(rocket, x_current, t_current, launchRampLength, initialDirection)
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
    if np.dot(position, initialDirection) <= launchRampLength + rocket.getHeight():
        stillAtLaunchRamp = True
    else:
        stillAtLaunchRamp = False
    # dPosition and dQuaternion
    RotationBody2Inertial = Kinematics.Rquaternion(q)
    RotationInertial2Body = RotationBody2Inertial.T
    dPosition = np.dot(RotationBody2Inertial,linearVelocity.T)
    dQuaternion = np.dot(Kinematics.quaternionGradient(q),angularVelocity.T)
    # forces in the body frame
    thrust = Forces.getThrust(rocket, t)
    gravity = np.dot(RotationInertial2Body, Forces.gravity(rocket,t).T)
    drag = rocket.SAMdrag(rocket, position, linearVelocity)
    AoA = 0
    lift = np.array([0,0,0])
    # inertia matrix and coriolis matrix for equations of motion
    # seen from origin of body frame, not from center of mass (See Fossen)
    H = Kinematics.TransformationMatrix(rocket.getCOM(t))
    m = rocket.getMass(t)
    I = rocket.getInertiaMatrix(t)
    IBody = splinalg.block_diag(m,m,m,I)
    IBody = np.dot(H.T, IBody)
    IBody = np.dot(IBody, H)
    S1 = Kinematics.CrossProductMatrix(m*angularVelocity)
    S2 = Kinematics.CrossProductMatrix(np.dot(I,angularVelocity.T))
    CBody = splinalg.block_diag(S1, -S2)
    CBody = np.dot(H.T,CBody)
    CBody = np.dot(CBody, H)
    # obtain generalized forces seen from origin of body frame
    totalForce = thrust + gravity + drag + lift
    if stillAtLaunchRamp:
        totalForce = np.dot(totalForce, initialDirection)*initialDirection
        totalMoment = np.array([0, 0, 0])
    else:
        arm = rocket.getCOP(AoA) - rocket.getCOM(t)
        totalMoment = np.cross(arm, drag + lift)
    genForceBody = np.concatenate((totalForce, totalMoment))
    genForceBody = np.dot(H.T, genForceBody)
    # find dx
    genVelocity = np.concatenate((linearVelocity, angularVelocity))
    rhs = genForceBody - np.dot(CBody, genVelocity.T)
    dGeneralizedVelocity = np.linalg.solve(IBody, rhs)
    if stillAtLaunchRamp: #defensive programming: this should be redundant
        dGenVelocity[4:6] = np.array([0, 0, 0])
        dQuaternion = np.array([0, 0, 0])
    dx = np.concatenate((dPosition, dQuaternion, dGenVelocity))
    return (dx, AoA, thrust, gravity, drag, lift)

def unwrapState(x):
    position = x[:,0:3]
    quaternion = x[:,3:7]
    linearVelocity = x[:,7:10]
    angularVelocity = x[:,10:13]
    return (position, quaternion, linearVelocity, angularVelocity)
