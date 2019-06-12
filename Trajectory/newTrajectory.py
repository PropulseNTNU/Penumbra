import sys
sys.path.append('../Forces/')
import Wind
import Kinematics
import Forces
from Rocket1 import Airbrakes
import numpy as np
import scipy.linalg as splinalg
from scipy.constants import g
import matplotlib.pyplot as plt

rad2deg = 180/np.pi
deg2rad = 1/rad2deg

def initialState(rocket, initialInclination):
    # Initial inclination of rocket
    initialPitch = np.pi/2-initialInclination
    # Calculate rotation matristate (the rocket orientation measured in world coords.)
    # at this initial inclination
    R = Kinematics.Ryzx(initialPitch, 0, 0)
    # rocket longitudinal astateis in world frame (state-astateis of body)
    initialDirection= R[:,0]
    initialPosition = rocket.getLength()*initialDirection
    initialQuaternion = Kinematics.euler2quaternion(initialPitch, 0, 0)
    initialLinearVelocity = np.array([0, 0, 0])
    initialAngularVelocity = np.array([0, 0, 0])
    # Initial state vector
    state0 = np.concatenate((initialPosition, initialQuaternion,
    initialLinearVelocity, initialAngularVelocity))
    return (state0, initialDirection)

def unwrapState(state):
    """
    A simple routine to unwrap the state vector to its costituent parts (position, orientation, velocity etc.)
    """
    position = state[0:3]
    quaternion = state[3:7]
    linearVelocity = state[7:10]
    angularVelocity = state[10:13]
    euler = Kinematics.quaternion2euler(quaternion)
    return (position, euler, linearVelocity, angularVelocity)

class Trajectory:
    def __init__(self, rocket, initialInclination, rampLength, simTime, dragDeviation = 0,
    windObj = Wind.nullWind(), air_brakes=Airbrakes(0,0,0)):
        # Integrator specs
        self.dt = 0.04 # temp
        self.timeArray = np.arange(0, simTime, self.dt)
        self.steps = len(self.timeArray)
        # Rocket object
        self.rocket = rocket
        self.airbrakes = air_brakes
        # Initial state
        self.initialInclination = initialInclination
        self.initState, self.initialDirection = initialState(rocket, initialInclination)
        self.rampLength = rampLength
        # State vector and Trajectory data
        self.state = self.initState
            # Translation
        self.position = np.zeros((self.steps, 3))
        self.linearVelocity = np.zeros((self.steps, 3))
        self.linearAcceleration = np.zeros((self.steps, 3))
            # Rotation
        self.orientation = np.zeros((self.steps, 3))
        self.AoA = np.zeros(self.steps)
        self.angularVelocity = np.zeros((self.steps, 3))
        self.angularAcceleration = np.zeros((self.steps, 3))
            # Forces
        self.drag = np.zeros((self.steps, 3))
        self.lift = np.zeros((self.steps, 3))
        self.gravity = np.zeros((self.steps, 3))
        self.thrust = np.zeros((self.steps, 3))
        # Environment
        self.wind = windObj
        self.dragDeviation = dragDeviation
        self.windVelocities = np.zeros((self.steps, 3))
        # Auxiliary
        self.simulationRun = False

    # Main functions
    def run(self, printStatistics = True):
        # Run Integrator
        print('Simulating..')
        self.__RK4()
        # Print message for status
        print('Simulation complete!')
        self.simulationRun = True
        if printStatistics:
            self.printTrajectoryStatistics()

    def flush(self):
        # flush trajectory data
        self.__init__(self.rocket, self.initialInclination, self.rampLength,
                      self.timeArray[-1]+self.dt, self.dragDeviation,
                      self.wind)
        self.simulationRun = False

    def plot(self):
        # Position
        t = self.timeArray
        position = self.position
        linearVelocity = self.linearVelocity
        angularVelocity = self.angularVelocity
        acceleration = self.linearAcceleration
        angularAcceleration = self.angularAcceleration
        orientation = self.orientation
        AoA = self.AoA
        thrust = self.thrust
        drag = self.drag
        lift = self.lift
        gravity = self.gravity
        ###
        plt.figure()
        ax1 = plt.subplot(311, ylabel='x [m]')
        ax1.plot(t, position[:,0], label='x', lw=2, c='r')
        ax1.grid()
        plt.subplots_adjust(hspace=0.5)
        ax2 = plt.subplot(312, ylabel='y [m]')
        ax2.plot(t, position[:,1], label='y', lw=2, c='b')
        ax2.grid()
        plt.subplots_adjust(hspace=0.5)
        ax3 = plt.subplot(313, xlabel='time [s]', ylabel=' altitude [m]')
        ax3.plot(t, -position[:,2], label='h', lw=2, c='g')
        ax3.grid()
        plt.figure()
        #ax1 = plt.subplot(111, xlabel='lat. distance [m]', ylabel=' altitude [m]', xlim=(0,1000), ylim=(0,3100))
        plt.plot(np.sqrt(position[:,0]**2+position[:,1]**2), -position[:,2], lw=3, c='r')
        plt.plot([1800],[0], lw=0)
        plt.title('altitude vs. lateral distance')
        plt.xlabel('lat. distance [m]')
        plt.ylabel('altitude [m]')
        plt.grid()
        plt.gcf().gca().set_aspect('equal')

        # Orientation
        plt.figure()
        plt.title('Orientation')
        ax1 = plt.subplot(411, ylabel='pitch [deg]')
        ax1.plot(t, orientation[:,0]*rad2deg, label='pitch', lw=2, c='r')
        ax1.grid()
        plt.subplots_adjust(hspace=0.5)
        ax2 = plt.subplot(412, ylabel='yaw [deg]')
        ax2.plot(t, orientation[:,1]*rad2deg, label='yaw', lw=2, c='b')
        ax2.grid()
        plt.subplots_adjust(hspace=0.5)
        ax3 = plt.subplot(413, ylabel='roll [deg]')
        ax3.plot(t, orientation[:,2]*rad2deg, label='roll', lw=2, c='g')
        ax3.grid()
        plt.subplots_adjust(hspace=0.5)
        ax1 = plt.subplot(414, xlabel='time [s]', ylabel='AoA [deg]')
        ax1.plot(t[15:], AoA[15:]*rad2deg, label='AoA', lw=2, c='k')
        ax1.grid()

        # Linear Velocity
        plt.figure()
        ax1 = plt.subplot(311, ylabel='v_x [m/s]')
        ax1.plot(t, linearVelocity[:,0], lw=2, c='r')
        ax1.grid()
        ax1.set_title('Velocity (world coords)')
        plt.subplots_adjust(hspace=0.5)
        ax2 = plt.subplot(312, ylabel='v_y [m/s]')
        ax2.plot(t, linearVelocity[:,1], lw=2, c='b')
        ax2.grid()
        plt.subplots_adjust(hspace=0.5)
        ax3 = plt.subplot(313, xlabel='time [s]', ylabel='v_z [m/s]')
        ax3.plot(t, -linearVelocity[:,2], lw=2, c='g')
        ax3.grid()

        # Forces
        plt.figure()
        ax1 = plt.subplot(311, ylabel='[N]')
        ax1.plot(t, thrust[:,0], label='thrust-x', lw=2, c='r')
        ax1.plot(t, drag[:,0], label='drag-x', lw=2, c='b')
        ax1.plot(t, lift[:,0], label='lift-x', lw=2, c='g')
        ax1.plot(t, gravity[:,0], label='gravity-x', lw=2, c='k')
        ax1.set_title('Forces (body-coords)')
        ax1.legend(loc='upper right')
        ax1.grid()
        plt.subplots_adjust(hspace=0.5)
        ax2 = plt.subplot(312, ylabel='[N]')
        ax2.plot(t, thrust[:,1], label='thrust-y', lw=2, c='r')
        ax2.plot(t, drag[:,1], label='drag-y', lw=2, c='b')
        ax2.plot(t, lift[:,1], label='lift-y', lw=2, c='g')
        ax2.plot(t, gravity[:,1], label='gravity-y', lw=2, c='k')
        ax2.grid()
        ax2.legend(loc='upper right')
        plt.subplots_adjust(hspace=0.5)
        ax3 = plt.subplot(313, xlabel='time [s]', ylabel='[N]')
        ax3.plot(t, thrust[:,2], label='thrust-z', lw=2, c='r')
        ax3.plot(t, drag[:,2], label='drag-z', lw=2, c='b')
        ax3.plot(t, lift[:,2], label='lift-z', lw=2, c='g')
        ax3.plot(t, gravity[:,2], label='gravity-z', lw=2, c='k')
        ax3.legend(loc='upper right')
        ax3.grid()

        #Acceleration
        plt.figure()
        ax1 = plt.subplot(311, ylabel='a_x [m/s^2]')
        ax1.plot(t, acceleration[:,0], lw=2, c='r')
        ax1.grid()
        ax1.set_title('Acceleration (world coords)')
        plt.subplots_adjust(hspace=0.5)
        ax2 = plt.subplot(312, ylabel='a_y [m/s^2]')
        ax2.plot(t, acceleration[:,1], lw=2, c='b')
        ax2.grid()
        plt.subplots_adjust(hspace=0.5)
        ax3 = plt.subplot(313, xlabel='time [s]', ylabel='a_z [m/s^2]')
        ax3.plot(t, -acceleration[:,2], lw=2, c='g')
        ax3.grid()

        # Angular velocity
        plt.figure()
        ax1 = plt.subplot(311, ylabel='w_x [rad/s]')
        ax1.plot(t, angularVelocity[:,0], label='roll rate', lw=2, c='r')
        ax1.grid()
        ax1.legend(loc='upper right')
        plt.subplots_adjust(hspace=0.5)
        ax2 = plt.subplot(312, ylabel='w_y [rad/s]')
        ax2.plot(t, angularVelocity[:,1], label='pitch rate', lw=2, c='b')
        ax2.grid()
        ax2.legend(loc='upper right')
        plt.subplots_adjust(hspace=0.5)
        ax3 = plt.subplot(313, xlabel='time [s]', ylabel='w_z [rad/s]')
        ax3.plot(t, angularVelocity[:,2], label='yaw rate', lw=2, c='g')
        ax3.grid()
        ax3.legend(loc='upper right')
        plt.show()

    def printTrajectoryStatistics(self):
        """
        :param rocket: [rocket object]
        :param position: [np.array] The position matrix in worldcoords (len(time)x3 matrix)
        :param velocity: [np.array] The velocity matrix in worldcoords (len(time)x3 matrix)
        :param time: [np.array] The time array
        """
        if not self.simulationRun:
            print('Trajectory data is empty! Please run a simulation and try again..')
            exit(1)

        x, y, z = self.position[:,0], self.position[:,1], self.position[:,2]
        vx, vy, vz = self.linearVelocity[:,0], self.linearVelocity[:,1], self.linearVelocity[:,2]

        # Statistics
        index_apogee = np.argmax(-z)
        index_max_vertical_speed = np.argmax(-vz)
        apogee = -z[index_apogee]
        apogee_time = self.timeArray[index_apogee]
        lat_distance = np.max(np.sqrt(x**2 + y**2))
        max_vertical_speed = -vz[index_max_vertical_speed]
        time_at_maxspeed = self.timeArray[index_max_vertical_speed]
        burn_out = self.rocket.getMotor().getBurnTime()
        time_1_10_deployed = -0.11 + self.airbrakes.Time_brakes()
        time_1_4_deployed  = -0.055 + self.airbrakes.Time_brakes()
        time_1_2_deployed = self.airbrakes.Time_brakes()
        v_1_10_deployed = np.linalg.norm(self.linearVelocity[int(time_1_10_deployed/self.dt)])
        v_1_4_deployed = np.linalg.norm(self.linearVelocity[int(time_1_4_deployed/self.dt)])
        v_1_2_deployed = np.linalg.norm(self.linearVelocity[int(time_1_2_deployed/self.dt)])

        print('\nTRAJECTORY STATISTICS:')
        print(33*'-')
        print('Max vertical speed: %1.1f m/s' %max_vertical_speed)
        print('\t Time at max speed: %1.1f s' %time_at_maxspeed)
        print('\t Time at burnout: %1.1f s' %burn_out)
        print('Air brakes info:')
        print('\tspeed at 10%% deployment: %1.1f m/s' %v_1_10_deployed)
        print('\tspeed at 25%% deployment: %1.1f m/s' %v_1_4_deployed)
        print('\tspeed at 50%% deployment: %1.1f m/s' %v_1_2_deployed)
        print('Apogee: %1.0f m' %apogee)
        print('\t Time at apogee: %1.1f s' %apogee_time)
        print('Lateral distance traveled: %1.0f m' %lat_distance)

    # Private functions here
    def __eqOfMotion(self, state, t):
        rocket = self.rocket
        epsilon = 1e-10
        position = state[0:3]
        quaternion = state[3:7]
        linearVelocity = state[7:10]
        angularVelocity = state[10:13]
        stillAtLaunchRamp = False
        # determine if whether at launch ramp or not
        if np.dot(position, self.initialDirection) <= self.rampLength + rocket.getLength():
            stillAtLaunchRamp = True
        else:
            stillAtLaunchRamp = False
        # dPosition and dQuaternion
        RotationBody2Inertial = Kinematics.Rquaternion(quaternion)
        RotationInertial2Body = RotationBody2Inertial.T
        # Velocity of rocket in world frame
        dPosition = RotationBody2Inertial @ linearVelocity
        dQuaternion = Kinematics.quaternionGradient(quaternion) @ angularVelocity
        # Add wind to current rocket velocity to get total air velocity
        windVelocity = self.wind.getWindVector(-position[2], t)
        airVelocity = windVelocity - dPosition
        airSpeed = np.linalg.norm(airVelocity)
        xAxisBody = RotationBody2Inertial[:,0]
        dirWindVelocity = (airVelocity/(airSpeed + epsilon))
        # definition of angle of attack
        AoA = np.arccos(np.dot(-dirWindVelocity, xAxisBody))
        # unit vector that points in drag direction (body coords.)
        dirDragBody = RotationInertial2Body @ (dirWindVelocity)
        projectedDragBody = np.array([0, dirDragBody[1], dirDragBody[2]])
        dirProjectedDragBody = projectedDragBody/(np.linalg.norm(projectedDragBody) + epsilon)
        # unit vector that points in lift direction (body coords.)
        dirLiftBody = np.sin(AoA)*np.array([1, 0, 0]) + np.cos(AoA)*dirProjectedDragBody
        # Forces in body coords
        aeroForces = rocket.getAeroForces(position, airVelocity, AoA)
        Cd_airbrakes = self.airbrakes.drag_coeff(position, t)
        drag = aeroForces[0]*dirDragBody + Cd_airbrakes*airSpeed**2*dirDragBody
        lift = aeroForces[1]*dirLiftBody
        thrust = np.array([rocket.getMotor().thrust(t), 0, 0])
        # gravity in body coords
        gravityWorld = np.array([0, 0, rocket.getMass(t)*g])
        gravityBody = RotationInertial2Body @ gravityWorld
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
        # If rocket is on launch ramp, don't allow it to rotate, only accelerate (fixed to ramp)
        if stillAtLaunchRamp:
            totalForce = np.array([totalForce[0], 0, 0])
            totalMoment = np.array([0, 0, 0])
        else:
            # After launch ramp, allow it to rotate (now calculating torques about COM)
            armAero = rocket.getCOP(position, airVelocity, AoA) - rocket.getCOM(t)
            totalMoment = np.cross(armAero, drag + lift)
        genForceBody = H.T @ np.concatenate((totalForce, totalMoment))
        # find dx
        genVelocity = np.concatenate((linearVelocity, angularVelocity))
        rhs = genForceBody - CBody @ genVelocity
        dGeneralizedVelocity = np.linalg.solve(IBody, rhs)
        dstate = np.concatenate((dPosition, dQuaternion, dGeneralizedVelocity))

        return dstate, AoA, forceMatrix, windVelocity

    def __RK4(self):
        """
        Runge-Kutta ODE solver of order 4, solving the equation system given by RHS
        :return: np.array[] ; matrix that contains the states at every instance (as row vectors)
        """
        # Initialize with initial conditions.
        state = self.initState

        # Runge-Kutta algorithm
        for i in range(self.steps):
            t = self.timeArray[i]
            s1, AoA, force, windVelocity = self.__eqOfMotion(state, t)
            s2 = self.__eqOfMotion(state + self.dt/2*s1, t + self.dt/2)[0] # get dx only
            s3 = self.__eqOfMotion(state + self.dt/2*s2, t + self.dt/2)[0] # get dx only
            s4 = self.__eqOfMotion(state + self.dt*s3, t + self.dt)[0] # get dx only

            dstate = (s1 + 2*s2 + 2*s3 + s4)/6
            state = state + self.dt*dstate
            self.__writeData(state, dstate, force, AoA, windVelocity, i)

    def __writeData(self, state, dstate, force, AoA, windvelocity, time_index):
        pos, orien, linvel, angvel = unwrapState(state)
        RotationBody2Inertial = Kinematics.Ryzx(orien[0], orien[1], orien[2])
        self.position[time_index] = pos
        self.linearVelocity[time_index] = RotationBody2Inertial @ linvel
        self.linearAcceleration[time_index] = RotationBody2Inertial @ dstate[7:10]
            # Rotation
        self.orientation[time_index] = orien
        self.AoA[time_index] = AoA
        self.angularVelocity[time_index] = angvel
        self.angularAcceleration[time_index] = dstate[10:13]
            # Forces
        self.drag[time_index] = force[:, 0]
        self.lift[time_index] = force[:, 1]
        self.gravity[time_index] = force[:, 2]
        self.thrust[time_index] = force[:, 3]
            # Wind
        self.windVelocities[time_index] = windvelocity

    # get functions here

    # set functions here
