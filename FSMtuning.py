"""
Script for running one complete simulation

Version: WIP
Last edit: 08.02.2019

--Propulse NTNU--
"""

import sys
sys.path.append('Rocket/')
sys.path.append('Forces/')
sys.path.append('Trajectory/')
sys.path.append('Visual')
import numpy as np
import scipy.linalg as splinalg
from Rocket1 import RocketSimple
import Trajectory
import Kinematics
import Forces
import pen_sensor as ps
import teensy_interface as ti
import matplotlib.pyplot as plt
import time
import FSMplotting as FSMplot

# Avoid division by 0 by adding epsilon to all denominators
epsilon = 1e-10

# Initialize a rocket
# FOR ROCKET CLASS 1
path1 = 'Tests/myRocket1/'
init_file1 = 'myRocket.dot'
Rocket = RocketSimple.from_file(init_file1, path1)
#Rocket.compressibleFlow(False)
Cd = 0.5
Across = np.pi*(Rocket.getBody().getDiameter()/2)**2

# Specify initial conditions
initialInclination = 6/180.0*np.pi
launchRampLength = 2.0*Rocket.getLength()
dt = 0.03
simulationTime= 30
x0, initialState = Trajectory.initialState(Rocket, initialInclination)

# Initialize serial port

# array of time, initialize arrays to store state vectors (one state for each row)
timelist = np.arange(0, simulationTime + dt, dt)
stateMatrix = np.zeros((len(timelist), len(x0)))

def RHS(x, t):
    return equationsMotion(x, t, Rocket, launchRampLength, initialState)

def plotData(teensyData, timeData, sumTime, ser):
    for key, val in teensyData.items():
        val[1].append(ti.readFloatData(ser, prefix=key, lines=100))
    timeData.append(sumTime)
    
   
def main():
    ser = ti.initSerial("/dev/ttyACM0", 9600, 1)#"/dev/ttyACM0"     /dev/cu.usbmodem4739891

    Aold = Across # Inital area
    # Initialize with initial conditions.
    x = x0
    stateMatrix[0] = x
    steps = len(timelist)

    Aabs = []
    xs = []
    dxs = []

    Aab = 0

    ## Data from teensy
    teensyData = {
        "t_h": ("Height[m]", []),
        "t_a": ("Acceleration[m/s^2]", []),
        "est_v": ("Estimated velocity[m/s]", []),
        "est_h": ("Estimated height[m]", []),
        "c_s": ("Airbrakes area[m^2]", []),
        "iter": ("Iteration time Penumbra[s]", []),
        "vel": ("Actual velocity[m/s]", [])
        }
    #for plotting and writing 
    timeData = []
    linearVelocity=[]
    position=[]

    FSMplot.init(teensyData, cols=2)

    iterationTime = 0
    simulatedTime = 0

    for i in range(1, steps):
        start = time.time()
        ser.flushInput()
        simulatedTime += dt

        Aab = ti.readFloatData(ser, prefix='c_s',lines= 100)

        plotData(teensyData, timeData, simulatedTime, ser)

        # Recieve data from serial port
        Anew = Across + Aab

        # Update rocket with new data
        Rocket.setCd(Cd*Anew/Aold)

        # Calculate equations of motion
        t = timelist[i]
        s1 = RHS(x, t)
        s2 = RHS(x + dt / 2 * s1, t + dt / 2)
        s3 = RHS(x + dt / 2 * s2, t + dt / 2)
        s4 = RHS(x + dt * s3, t + dt)

        # Update state
        dx = (s1 + 2 * s2 + 2 * s3 + s4)/6
        x = x + dt*dx
        stateMatrix[i] = x  # Store new state
        Aold = Anew # Update old area for next iteration
        linearVelocity.append(x[7])
        position.append(x[2])

        # Acceleration
        quaternion = x[3:7]
        Rbody2Inertial = Kinematics.Rquaternion(quaternion)
        accWorld = Rbody2Inertial @ dx[7:10].T  

        #send data to teensy
        ti.sendData(ser, [-x[2], -accWorld[2], iterationTime/i, x[7]])

        Aabs = Aabs + [[Aab]]
        xs = xs + [[x[2]]]
        dxs = dxs + [[np.linalg.norm(accWorld)]]
        iterationTime += time.time() - start
    print("Average iteration time: ", iterationTime/steps)
    FSMplot.plotData(teensyData, timeData)
    plt.show()
    plt.pause(60*30)
    

    lookUpTable=[]
    hoydeN=0;
    hoyde=0;
    diff=0;
    for i in range(len(position)):
       hoydeN=-int(np.floor(position[i]))
       if hoydeN < hoyde:
           continue
       if hoydeN>hoyde:
           #var2=int(np.floor(position[:,2][i+1]))
           diff=hoydeN-hoyde
           if i==0:
               a=(linearVelocity[i])/diff
               for j in range(diff):
                   lookUpTable.append(a*j)
                   hoyde+=1
           else:
               a=(linearVelocity[i]-(linearVelocity[i-1]))/diff
               for j in range(diff):
                   lookUpTable.append(a*j+lookUpTable[hoyde-j-1])
                   hoyde+=1
       lookUpTable.append(linearVelocity[i])
       hoyde+=1
    print(lookUpTable)



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
    windVelocity = np.array([0, 0, 0])
    # Add wind to current rocket velocity to get total air velocity
    airVelocity = dPosition + windVelocity
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

    return dx


if __name__ == '__main__':
    main()
