"""
The definition of the CFD rocket

Version: WIP
Last edit: 12.06.2019

--Propulse NTNU--
"""
import sys
sys.path.append('../Forces/')
import Forces as Forces
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp2d
from Rocket1 import Motor
from lib.File_utilities import find_parameter, unwrap_CFD_report

# Define some things for plotting
font = {'family': 'sans-serif', 'weight': 'bold', 'size': 16}
plt.rc('font', **font)
plt.rcParams['text.latex.preamble'] = [r'\boldmath']

class RocketCFD:
    def __init__(self, *args):
        print('Initializing rocket..')
        self.__initMass = args[0]
        self.__initInertiaMatrix = args[1]
        self.__initCOM = args[2]
        self.__length = args[3]
        self.__motor = args[4]

        # Print motor specs
        print(self.__motor)

        # for plotting, array of time during burn phase
        self.__time = np.arange(0, self.__motor.getBurnTime(), 5e-4)

        # Interpolation of forces and moments
        self.__drag = args[5]
        self.__lift = args[6]
        self.__moment = args[7]
        self.__freeAirStreamSpeeds = args[8]
        self.__AoAarray = args[9]
# =============================================================================
#         # Create rotation matrix for all AoA
#         def rotationMatrix(AoA):
#             sa, ca = np.sin(AoA*np.pi/180.0), np.cos(AoA*np.pi/180.0)
#             return np.array([[ca, sa],[-sa, ca]])
#      
#         # Initialize aeroforces (a 2xn matrix with drag as row 1 and lift as row 2)
#         for i in range(len(self.__AoAarray)):
#             rotMat = rotationMatrix(self.__AoAarray[i])
#             self.__aeroForces[:,i] = rotMat@self.__aeroForces[:,i]
# =============================================================================
        print('\tInterpolating the drag force..')
        self.__Dragforce = interp2d(self.__AoAarray, self.__freeAirStreamSpeeds, self.__drag)
        print('\tInterpolating the lift force..')
        self.__Liftforce = interp2d(self.__AoAarray, self.__freeAirStreamSpeeds, self.__lift)
        print('\tInterpolating the moment about COM (component normal to aerodynamic plane)..')
        self.__MomentAboutCOM = interp2d(self.__AoAarray, self.__freeAirStreamSpeeds, self.__moment)
        # Done
        print('Rocket initialized!\n')

    def getMotor(self):
        return self.__motor

    def getMass(self, t):
        return self.__initMass + self.__motor.getMass(t) - self.__motor.getMass(0)

    def getCOM(self, t):
        """
        :param t: [float] at time t [sec]

        :return: [float] the position of COM relative to nose tip [m]
        """
        mass = self.getMass(t)
        initMass = self.__initMass
        initCOM = self.__initCOM
        rocketLength = self.__length
        motorHeight = self.getMotor().getLength()
        initMotorMass = self.getMotor().getMass(0)
        motorMass = self.getMotor().getMass(t)
        COM = 1/mass*(initMass*initCOM + (rocketLength - motorHeight/2)*(initMotorMass - motorMass))

        return np.array([COM, 0, 0])

    def getLength(self):
        """
        :return: the length of the rocket [m]
        """
        return self.__length

    # Aero dynamics
    def getAeroForces(self, position, velocity, AoA):
        """
        :param velocity: [np.array] the air velocity relative to rocket [m/s]
        :param position: [np.array] The position vector in world coordinates
        :param AoA: [float] the angle of attack [rad]

        :return: [np.array] ([drag, lift]) on rocket attacking in COP [N]
        """
        z = abs(position[2])  # Vertical position of rocket
        speed = np.linalg.norm(velocity)
        density_reduction = np.exp(-z/Forces.h)  # Account for decreasing air density
        drag = self.__Dragforce(AoA, speed/Forces.c)
        lift = self.__Liftforce(AoA, speed/Forces.c)
        return np.array([drag, lift])*density_reduction

    def getMomentAboutCOM(self, position, velocity, AoA):
        """
        :param speed: [float] the air speed relative to rocket [m/s]
        :param position: [np.array] The position vector in world coordinates
        :param AoA: [float] the angle of attack [rad]

        :return: [float] The total moment on rocket about COM (component normal to aerodynamic plane) [Nm]
        """
        z = abs(position[2])  # Vertical position of rocket
        speed = np.linalg.norm(velocity)
        density_reduction = np.exp(-z/Forces.h)  # Account for decreasing air density
        return self.__MomentAboutCOM(AoA, speed)*density_reduction

    def getCOP(self, position, velocity, AoA):
        """
        :param position: [np.array] The position vector in world coordinates
        :param AoA: [float] the angle of attack [rad]

        :return: [float] The position of COP relative to nose tip [m]
        """
        F = self.getAeroForces(position, velocity, AoA)
        Fd = F[0]
        Fl = F[1]
        M = self.getMomentAboutCOM(position, velocity, AoA)
        COM_0 = self.getCOM(0)
        COPx = COM_0[0] - M/(Fl*np.cos(AoA) + Fd*np.sin(AoA) + 1e-10)
        COPx = COPx[0]
        return np.array([COPx, 0, 0])

    def getAeroData(self):
        return self.__AoAarray, self.__freeAirStreamSpeeds, self.__drag, self.__lift, self.__moment
    
    def getStabilityMargin(self, position, velocity, AoA, t=0):
        """
        Stability margin in units of body diameters (body calibers)
        """
        COM = self.getCOM(t)[0]
        COP = self.getCOP(position, velocity, AoA)[0]
        return COM - COP

    def compressibleFlow(self, state):
        if type(state) == bool:
            self.__compressibility = state
        else:
            print("Error: Enter either True or False in 'compressibleFlow'.")
            exit(1)

    def getCompressibilityState(self):
        return self.__compressibility
    
    def getInertiaMatrix(self, t):
        """
        :param t: [float] point in time [s]
        :return [np.array, 3x3 matrix] The inertia matrix at time t, centered at body origin [kgm^2]

        NOTE: Currently assuming that fuel is burning radially.
        """
        # TODO Add inertia for axial burning fuel.
        I0 = self.__initInertiaMatrix
        rInitMass = self.__initMass
        rInitCOMx = self.__initCOM
        rCOMx = self.getCOM(t)[0]
        mInitMass = self.__motor.getMass(0)
        mInitCOMx = self.__motor.getCOM(0)[0]
        mInitInertia = self.__motor.getInertiaMatrix(0)
        mMass = self.__motor.getMass(t)
        mInertia = self.__motor.getInertiaMatrix(t)
        deltaR = rCOMx - rInitCOMx
        deltaM = rCOMx + mInitCOMx + self.getLength()
        m = self.getMass(t) - mMass
        return I0 + rInitMass*np.diag([0, deltaR**2, deltaR**2]) + (mInertia - mInitInertia) + (
                    mMass - mInitMass)*np.diag([0, deltaM**2, deltaM**2]) + m*rCOMx**2*np.diag([0, 1, 1])

    # auxiliary
    def plot(self):
        # Plots of motor performance
        self.__motor.plotPerformance(False)  # plt.show=False in this case

        # Plots of rocket characteristics (during burn phase)
        time = self.__time
        burnTime = self.__motor.getBurnTime()
        print('Creating plot of COM and mass...')
        COM = np.array([self.getCOM(t)[0] for t in time])
        mass = np.array([self.getMass(t) for t in time])
        plt.figure()
        ax1 = plt.subplot(211, xlabel='time [s]', ylabel='[cm]')
        ax1.plot(time, 100*COM, label='COMx(t)', lw=2, c='r')
        ax1.set_title('COM of rocket during burn phase of %1.1f s'%burnTime)
        ax1.grid()

        ax2 = plt.subplot(212, xlabel='time [s]', ylabel='[g]')
        ax2.plot(time, mass*1e3, label='mass(t)', lw=2, c='b')
        ax2.set_title('Mass of rocket')
        ax2.grid()
        plt.subplots_adjust(hspace=0.5)
        # Moment of inertia
        print('Creating plots of rotational/longitudinal moment of inertia...')
        steps = len(time)
        Ixx = np.zeros(steps)
        Iyy = np.zeros(steps)
        for i in range(steps):
            MOI = np.diagonal(self.getInertiaMatrix(time[i]))
            Ixx[i] = MOI[0]
            Iyy[i] = MOI[1]
        # Izz = Iyy with principle axes, don't need to create another array for this
        plt.figure()
        ax1 = plt.subplot(211, xlabel='time [s]', ylabel='[kgcm²]')
        ax1.plot(time, Ixx*1e4, label='Ixx(t)', lw=2, c='r')  # Multiplied by 10,000 to get correct unit
        ax1.set_title('Rotational inertia of rocket during burn phase of %1.1f s'%burnTime)
        ax1.grid()

        ax2 = plt.subplot(212, xlabel='time [s]', ylabel='[kgcm²]')
        ax2.plot(time, Iyy*1e4, label='I(t)', lw=2, c='b')  # Multiplied by 10,000 to get correct unit
        ax2.set_title('Longitudinal inertia')
        ax2.grid()
        plt.subplots_adjust(hspace=0.5)
        plt.show()
        print('Plotting done!\n')

    
    @staticmethod
    def from_file(initFile, drag_report, lift_report, moment_report, path_to_file=''):
        """
                Create an instance of CFDrocket by reading some files

                example: initFile.dot as initFile, myMotor.dot as motor and CFD files as drag.dot, lift.dot, moments.dot

                        in 'initFile.dot':

                        initial_mass = value of initial mass
                        initial_moi = Ixx, Iyy, Izz
                        initial_com = x-value of initial COM (from nosetip)
                        length = value of rocket total length
                         motor = myMotor.dot

                        myRocket = Rocket.from_file('initFile.dot', 'drag_report', 
                                                    'lift_report', 'moment_report')

                :param path_to_file: The path to the files above relative to current path (none by default)

                return: A rocket instance with specs from initFile and CFD.
                """
        path = path_to_file + initFile
        initMass = find_parameter(path, 'initial_mass')  # in kilo grams
        initMOI = find_parameter(path, 'initial_moi')
        initMOI = np.diag(np.array([x.strip() for x in initMOI.split(',')]).astype(float))  # in kg*m^2
        initCOM = find_parameter(path, 'initial_com')  #in meters
        length = find_parameter(path, 'length')  #in meters
        motor = Motor.from_file(path_to_file + find_parameter(path, 'motor'))

        path1 = path_to_file + drag_report
        path2 = path_to_file + lift_report
        path3 = path_to_file + moment_report

        drag, speed, AoA = unwrap_CFD_report(path1)
        lift = unwrap_CFD_report(path2)[0]
        moment = unwrap_CFD_report(path3)[0]

        return RocketCFD(float(initMass), initMOI, float(initCOM), float(length), motor, drag, lift, moment, speed, AoA)
