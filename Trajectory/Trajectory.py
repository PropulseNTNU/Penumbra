import numpy as np
import scipy.linalg as splinalg
import scipy.integrate as spintegrate
import Rocket.Rocket as RockMod
import Force.Forces as Forces

# rotations, kinematics, kinetics help functions
def Ryzx(pitch, yaw, roll):
    cp = np.cos(pitch)
    sp = np.sin(pitch)
    cy = np.cos(yaw)
    sy = np.sin(yaw)
    cr = np.cos(roll)
    sr = np.sin(roll)
    R = np.array([[cp*cy, sp*sr-cp*sy*cr, sp*cr+cp*sy*sr],
    [sy, cy*cr, -cy*sr], [-sp*cy, cp*sr+sp*sy*cr, cp*cr-sp*sy*sr]])
    return R

def Rquaternion(q):
    q = q/np.linalg.norm(q)
    eta = q[0]
    eps = q[1:4]
    S = CrossProductMatrix(eps)
    R = np.eye(3) + 2*eta*S + 2*np.dot(S,S)
    return R

def euler2quaternion(pitch, yaw, roll):
    R = Ryzx(pitch, yaw, roll)
    temp = np.array([R[0,0], R[1,1], R[2,2], R.trace()])
    index = temp.argmax()
    print(index)
    p_i = np.sqrt(1 + 2*temp[index] - R.trace())
    if index==0:
        p1 = p_i
        p2 = (R[1,0]+R[0,1])/p_i
        p3 = (R[0,2]+R[2,0])/p_i
        p4 = (R[2,1]-R[1,2])/p_i
    elif index==1:
        p1 = (R[1,0]+R[0,1])/p_i
        p2 = p_i
        p3 = (R[2,1]+R[1,2])/p_i
        p4 = (R[0,2]-R[2,0])/p_i
    elif index==2:
        p1 = (R[0,2]+R[2,0])/p_i
        p2 = (R[2,1]+R[1,2])/p_i
        p3 = p_i
        p4 = (R[1,0]-R[0,1])/p_i
    else:
        p1 = (R[2,1]-R[1,2])/p_i
        p2 = (R[0,2]-R[2,0])/p_i
        p3 = (R[1,0]-R[0,1])/p_i
        p4 = p_i
    q = 0.5*np.array([p4, p1, p2, p3])
    q = q/np.linalg.norm(q)
    return q

def quaternion2euler(q):
    R = Rquaternion(q)
    pitch = np.atan2(-R[2,0], R[0,0])
    yaw = np.asin(R[1,0])
    roll = np.atan2(-R[1,2], R[1,1])
    return np.array([pitch, yaw, roll])

def quaternionGradient(q):
    q = q/np.linalg.norm(q)
    T = 0.5* np.array([[-q[1], -q[2], -q[3]], [q[0], -q[3], q[2]], [q[3], q[0], -q[1]], [-q[2], q[1], q[0]]])
    return T

def CrossProductMatrix(v):
    S = np.array([[0, -v[2], v[1]],[v[2], 0, -v[0]],[-v[1], v[0], 0]])
    return S

def TransformationMatrix(positionVector):
    S = CrossProductMatrix(positionVector)
    H = np.eye(6)
    H[0:3,3:6] = S
    return H

# trajectory module
def calculateTrajectory(rocket, initialInclination, launchRampLength, timeStep, simulationTime):
    # x is the state of the vector
    # x = [position, quaternion, linear velocity, angular velocity]
    (x0, initialDirection) = initialState(rocket, initialInclination)
    (t, x, AoA, thrust, gravity, drag, lift) = integrateEquationsMotion(rocket, x0, launchRampLength, initialDirection, timeStep, simulationTime)
    (position, quaternion, linearVelocity, angularVelocity) = unwrapState(x)
    return ((t, position, quaternion, linearVelocity, angularVelocity, AoA, thrust, gravity, drag, lift))

def initialState(rocket, initialInclination):
    R = Ryzx(initialInclination-90, 0, 0)
    initialDirection= R[:,0]
    initialPosition = rocket.getLength()*initialDirection
    initialQuaternion = euler2quaternion(initialInclination,0,0)
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

def equationsMotion(rocket, x, t, launchRampLength, initialDirection)
    position = x[0:3]
    quaternion = x[3:7]
    linearVelocity = x[7:10]
    angularVelocity = x[10:13]
    # determine if whether at launch ramp or not
    if np.dot(position, initialDirection)<=launchRampLength:
        stillAtLaunchRamp = True
    else:
        stillAtLaunchRamp = False
    # dPosition and dQuaternion
    RotationBody2Inertial = Rquaternion(q)
    RotationInertial2Body = RotationBody2Inertial.T
    dPosition = np.dot(RotationBody2Inertial,linearVelocity.T)
    dQuaternion = np.dot(quaternionGradient(q),angularVelocity.T)
    # forces in the body frame
    thrust = Forces.getThrust(rocket, t)
    gravity = np.dot(RotationInertial2Body, Forces.gravity(rocket,t).T)
    drag = rocket.SAMdrag(rocket, position, linearVelocity)
    AoA = 0
    lift = np.array([0,0,0])
    # inertia matrix and coriolis matrix for equations of motion
    # seen from origin of body frame, not from center of mass (See Fossen)
    H = TransformationMatrix(rocket.getCOM(t))
    m = rocket.getMass(t)
    I = rocket.getInertiaMatrix(t)
    IBody = splinalg.block_diag(m,m,m,I)
    IBody = np.dot(H.T, IBody)
    IBody = np.dot(IBody, H)
    S1 = CrossProductMatrix(m*angularVelocity)
    S2 = CrossProductMatrix(np.dot(I,angularVelocity.T))
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

def test():
    deg2rad = np.pi/180
    rad2deg = 180/np.pi
    R0 = Ryzx(0,0,0)
    R1 = Ryzx(np.pi/4,0,0)
    R2 = Ryzx(0,np.pi/4,0)
    R3 = Ryzx(0,0,np.pi/4)
    print(R0)
    print(R1)
    print(R2)
    print(R3)
    R4 = Ryzx(-80*deg2rad,0,0)
    print(R4)
    w = np.array([1.5, -2.3, 3.1])
    print(CrossProductMatrix(w))
    print(np.dot(CrossProductMatrix(w),w.T))
    myRocket = SimpleRocketModel()
    H = TransformationMatrix(myRocket.getCenterMass(0))
    print(H)
    print(IBody(myRocket,1))
    print(CoriolisMatrix(myRocket,2,np.array([0,0,0])))
    print(CoriolisMatrix(myRocket,2,np.array([2,-2,3])))
    print(rhsDynamicsBodyFrame(myRocket, 1, np.array([4,3,2,1,2,3])))
    pitch = np.pi/3
    yaw = -np.pi/7
    roll = np.pi/4
    R5 = Ryzx(pitch, yaw, roll)
    q = euler2quaternion(pitch, yaw, roll)
    R6 = Rquaternion(q)
    print(R5)
    print(R6)
    print(R5-R6)
    x=np.array([1., 2., 3., 4., 5., 6., 7., 8., 9., 10., 11., 12., 13.])
    t=2
    print(rhsEquationsOfMotion(x, t, myRocket))

def main():
    #test()
    Nose = Nose.from_file('test.dot', 'conic')
    Fin = Fin.from_file('fin.dot')
    Body = Body.from_file('body.dot')
    Motor = Motor.from_file('motor.dot')
    Payload = Payload.from_file('payload.dot')
    partsPlacement = np.array() #TODO

    rocket = Rocket.from_file('rocket.dot')

    initialInclination = 6.0/180.0*np.pi
    x0 = initialState(rocket, initialInclination)
    #trajectory = spintegrate.ode(rhsEquationsOfMotion).set_integrator('zvode', method='bdf')
    #trajectory.set_initial_value(x0,0).set_f_params(rocket)
    #tEnd = 200
    #dt = 0.01
    #while trajectory.successful() and trajectory.t < tEnd:
    #    trajectory.integrate(trajectory.t+dt)
    #    print("%g %g" % (trajectory.t, trajectory.y))
    t = np.linspace(0,100,1000)
    rhsODE = lambda x, t: rhsEquationsOfMotion(x, t, rocket)
    #print(rhsODE(x0, 0))
    x = spintegrate.odeint(rhsODE, x0, t)
    #x = spintegrate.odeint(rhsEquationsOfMotion, x0, t, args=(rocket))
    print(x)
main()
