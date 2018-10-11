import numpy as np
import scipy.linalg as splinalg
import scipy.integrate as spintegrate

g = 9.81
#Rocket module
class SimpleRocketModel:
    def __init__(self):
        self.height = 5
        self.diameter = 0.5
        self.dragCoefficient = 1
        self.liftCoefficient = 1

    def getTotalHeight(self):
        return self.height

    def getCrossSection(self):
        return self.diameter

    def getMass(self,t):
        return 1000

    def getCenterMass(self,t):
        return np.array([-3, 0, 0])

    def getInertiaMatrix(self, t):
        return np.array([[1,0,0], [0,1,0], [0,0,10]])

    def getThrust(self, t):
        return 1000*t*np.exp(-5*np.power(t,2))

    def getExhaustVelocity(self,t):
        return 1

    def getDragCoefficient(self, angleOfAttack):
        return self.dragCoefficient*np.sin(angleOfAttack)

    def getLiftCoefficient(self, angleOfAttack):
        return self.liftCoefficient*np.sin(angleOfAttack)

    def getCenterPressure(self, angleOfAttack):
        return np.array([-3, 0, 0])

# kinematics module
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

def quaternionGradientMatrix(q):
    q = q/np.linalg.norm(q)
    T = 0.5* np.array([[-q[1], -q[2], -q[3]], [q[0], -q[3], q[2]], [q[3], q[0], -q[1]], [-q[2], q[1], q[0]]])
    return T

def CrossProductMatrix(v):
    S = np.array([[0, -v[2], v[1]],[v[2], 0, -v[0]],[-v[1], v[0], 0]])
    return S

# kinetics module
def TransformationMatrix(positionVector):
    S = CrossProductMatrix(positionVector)
    H = np.eye(6)
    H[0:3,3:6] = S
    return H

def InertiaMatrix(rocket, t):
    H = TransformationMatrix(rocket.getCenterMass(t))
    m = rocket.getMass(t)
    I = rocket.getInertiaMatrix(t)
    InertiaMatrix = splinalg.block_diag(m,m,m,I)
    InertiaMatrix = np.dot(H.T, InertiaMatrix)
    InertiaMatrix = np.dot(InertiaMatrix, H)
    return InertiaMatrix

def CoriolisMatrix(rocket, t, angularVelocity):
    H = TransformationMatrix(rocket.getCenterMass(t))
    m = rocket.getMass(t)
    I = rocket.getInertiaMatrix(t)
    S1 = CrossProductMatrix(m*angularVelocity)
    S2 = CrossProductMatrix(np.dot(I,angularVelocity.T))
    CoriolisMatrix = splinalg.block_diag(S1, -S2)
    CoriolisMatrix = np.dot(H.T,CoriolisMatrix)
    CoriolisMatrix = np.dot(CoriolisMatrix, H)
    return CoriolisMatrix

def rhsDynamicsBodyFrame(rocket, t, velocities):
    angularVelocity = velocities[3:6]
    aerodynamicForce = np.array([0,0,0]) #TODO
    generalizedForcesOnCenterMass = np.array([rocket.getThrust(t),0,0,0,0,0])
    H = TransformationMatrix(rocket.getCenterMass(t))
    I = InertiaMatrix(rocket, t)
    C = CoriolisMatrix(rocket, t, angularVelocity)
    rhs = np.dot(H.T, generalizedForcesOnCenterMass)- np.dot(C, generalizedVelocity.T)
    dGeneralizedVelocity = np.linalg.solve(I, rhs)
    return dGeneralizedVelocity

def thrustForce(rocket, t):
    thrustForce = np.array([rocket.getThrust(t)*rocket.getExhaustVelocity(t), 0, 0])
    return thrustForce

def gravityForce(rocket, t, quaternion):
    RotationInertial2Body = Rquaternion(quaternion).T
    gravityForce = np.dot(RotationInertial2Body,np.array([0,0,-rocket.getMass(t)*g]))
    return gravityForce

def aerodynamicForceAndMoment(rocket, t, quaternion, linearVelocity, airDensity):
    RotationBody2Inertial = Rquaternion(q)
    windVelocity = np.array([0, 0, 0])
    airVelocity = linearVelocity - windVelocity
    xAxisBody = np.dot(RotationBody2Inertial,np.array([1, 0, 0]).T)
    angleOfAttack = np.arccos(np.dot(x_b, airVelocity)/np.linalg.norm(relativeVelocity))
    dragCoefficient = rocket.getDragCoefficient(angleOfAttack)
    liftCoefficient = rocket.getLiftCoefficient(angleOfAttack)
    commonFactor = 0.5*airDensity*np.power(airVelocity,2)*rocket.getCrossSection()
    #TODO

def funEquationsOfMotion(x, t, rocket):
    # Also Jacobian of rhs of equations of motion
    position = x[0:3]
    quaternion = x[3:7]
    RotationBody2Inertial = Rquaternion(q)
    linearVelocity = x[7:10]
    angularVelocity = x[10:13]
    dPosition = np.dot(RotationBody2Inertial,linearVelocity.T)
    dQuaternion = np.dot(quaternionGradientMatrix(q),angularVelocity.T)
    #TODO
    dragForce = rocket.getDragCoefficient
    liftForce =
    totalForce = gravityForce + thrustForce + dragForce + liftForce
    gravityMoment = np.array([0, 0, 0])
    thrustMoment = np.array([0, 0, 0])
    dragMoment =
    liftMoment =
    totalMoment = gravityMoment + thrustMoment + dragMoment + liftMoment
    generalizedForces = np.array([]) concatenate


    dv = rhsDynamicsBodyFrame(rocket, t, x[7:13])
    dx = np.concatenate((dp, dq, dv))
    return dx

def initialState(rocket, initialInclination):
    R = Ryzx(initialInclination-90, 0, 0)
    initialPosition = rocket.getTotalHeight()*R[:,0]
    initialQuaternion = euler2quaternion(initialInclination,0,0)
    initialLinearVelocity = np.array([0, 0, 0])
    initialAngularVelocity = np.array([0, 0, 0])
    x0 = np.concatenate((initialPosition, initialQuaternion,
    initialLinearVelocity, initialAngularVelocity))
    return x0

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
    print(InertiaMatrix(myRocket,1))
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
    rocket = SimpleRocketModel()
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
