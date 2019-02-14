"""
Rocket module - The definition of the rocket and its constituent parts

Version: WIP
Last edit: 04.02.2019

--Propulse NTNU--
"""
import sys
sys.path.append('../Forces/')
import Forces as Forces

import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from lib.File_utilities import find_parameter

# Define some things for plotting
font = {'family': 'sans-serif', 'weight': 'bold', 'size': 16}
plt.rc('font', **font)
plt.rcParams['text.latex.preamble'] = [r'\boldmath']

# Geometry type for nose
noseTypes = ['cone', 'hemisphere', 'ogive']

class Nose:
    # Constructor for conic nose
    def __init__(self, noseType, *args):
        # Member variables
        self.__noseType = noseType
        if noseType == noseTypes[0] or noseType == noseTypes[2]:  # Conic or Ogive
            self.__diameter = args[0]
            self.__length = args[1]
            self.__thickness = args[2]
            self.__density = args[3]
        elif noseType == noseTypes[1]:  # Hemisphere
            self.__diameter = args[0]
            self.__length = self.__diameter/2
            self.__thickness = args[1]
            self.__density = args[2]
        print("Nose initialized!\n")

    def __str__(self):
        D = str(self.__diameter)
        d = str(self.__thickness*1e3)
        m = str(round(self.getMass(), 2))
        rho = str(self.__density)
        if self.__noseType == noseTypes[0]:  # Conic
            h = str(self.__length)
            return "Diameter: " + D + " m\n" + "Height: " + h + " m\n" + "Thickness: " + d + " mm\n" + "Mass: " + m + \
                   " kg\n" + "Density: " + rho + " kgm^-3"
        elif self.__noseType == noseTypes[1]:  # Hemisphere
            return "Diameter: " + D + " m\n" + "Thickness: " + d + " mm\n" + "Mass: " + m + \
                   " kg\n" + "Density: " + rho + " kgm^-3"

    # Member functions
    def getSurfaceArea(self):
        if self.__noseType == noseTypes[0]:  # Conic
            R = self.__diameter/2  # Outer radius of cone
            H = self.__length
            L = np.sqrt(R**2 + H**2) # distance from cone edge to cone tip.
            return np.pi*R*L
        elif self.__noseType == noseTypes[1]:  # Hemisphere
            R = self.__diameter/2
            return 2*np.pi*(R**2)
        elif self.__noseType == noseTypes[2]:  # Ogive
            R = self.__diameter/2
            L = self.__length
            a = (R**2 + L**2)/(2*R) # Ogive parameter
            return 2*np.pi*((R-a)*np.arcsin(L/a) + L)

    def getVolume(self):
        if self.__noseType == noseTypes[0]:  # Conic
            d = self.__thickness
            R2 = self.__diameter/2  # Outer radius of cone
            H2 = self.__length
            H1 = H2 - d
            R1 = R2*H1/H2
            return np.pi/3*(H2*R2**2 - H1*R1**2)
        elif self.__noseType == noseTypes[1]:  # Hemisphere
            d = self.__thickness
            R2 = self.__diameter/2
            R1 = R2 - d
            return 2/3*np.pi*(R2**3 - R1**3)
        elif self.__noseType == noseTypes[2]:  # Ogive
            d = self.__thickness
            R2 = self.__diameter/2
            H2 = self.__length
            H1 = H2 - d
            R1 = R2*H1/H2
            return ((3*np.pi**2)/32)*R2**2*H2-((3*np.pi**2)/32)*R1**2*H1

    def getCavityVolum(self):
        if self.__noseType == noseTypes[0]:  # Conic
            d = self.__thickness
            R2 = self.__diameter/2  # Outer radius of cone
            H2 = self.__length
            H1 = H2 - d
            R1 = R2*H1/H2
            return np.pi/3*H1*R1**2
        elif self.__noseType == noseTypes[1]:  # Hemisphere
            d = self.__thickness
            R2 = self.__diameter/2
            R1 = R2 - d
            return 2/3*np.pi*R1**3
        elif self.__noseType == noseTypes[2]:  # Ogive
            d = self.__thickness
            R2 = self.__diameter/2
            H2 = self.__length
            H1 = H2 - d
            R1 = R2*H1/H2
            return ((3*np.pi**2)/32)*R1**2*H1

    def getLength(self):
        if self.__noseType == noseTypes[0] or self.__noseType == noseTypes[2]:
            return self.__length
        elif self.__noseType == noseTypes[1]:
            return self.__diameter/2

    def getDiameter(self):
        return self.__diameter

    def getMass(self):
        return self.__density*self.getVolume()

    def getInertiaMatrix(self):
        if self.__noseType == noseTypes[0] or self.__noseType == noseTypes[2]:  # Conic and Ogive
            r = self.__diameter/2
            m = self.getMass()
            Ixx = 1/2*m*r**2
            Iyy = 3*Ixx - m*self.getCOM()[0]**2
            Izz = Iyy
            return np.diag([Ixx, Iyy, Izz])
        elif self.__noseType == noseTypes[1]:  # Hemisphere
            r = self.__diameter/2
            m = self.getMass()
            Ixx = 2/3*m*r**2
            Iyy = 3*Ixx - m*self.getCOM()[0]**2
            Izz = Iyy
            return np.diag([Ixx, Iyy, Izz])

    def getNoseType(self):
        return self.__noseType

    def getCOM(self):
        COM = 0
        if self.__noseType == noseTypes[0]:
            V = self.getVolume()
            d = self.__thickness
            R2 = self.__diameter/2  # Outer radius of cone
            H2 = self.__length
            H1 = H2 - d
            R1 = R2*H1/H2
            a = 1/2*H2**2 - 2*H2**3/(3*H1) + H2**4/(4*H1**2)
            COM = np.pi/V*((H2*R2)**2/12 - a*R1**2)
        elif self.__noseType == noseTypes[1]:
            d = self.__thickness
            R2 = self.__diameter/2
            R1 = R2 - d
            COM = 3/8*(R2**4 - R1**4)/(R2**3 - R1**3)
        elif self.__noseType == noseTypes[2]:
            V = self.getVolume()
            d = self.__thickness
            R2 = self.__diameter/2  # Outer radius of cone
            H2 = self.__length
            H1 = H2 - d
            R1 = R2*H1/H2
            COM = ((32/(15*np.pi)) * (H2**2 * R2**2 - H1**2 * R1**2)
            /(H2 * R2**2 - H1 * R1**2) - (H2/2))
        # COM relative to bottom surface of nose
        return np.array([COM, 0, 0])

    @staticmethod
    def from_file(file):
        noseType = find_parameter(file, "nose_type")
        if noseType.lower() == noseTypes[0] or noseType.lower() == noseTypes[2]:  # Conic or Ogive
            diameter = find_parameter(file, "diameter")
            length = find_parameter(file, "length")
            density = find_parameter(file, "density")
            thickness = find_parameter(file, "thickness")
            return Nose(noseType.lower(), eval(diameter), eval(length), eval(thickness), eval(density))
        elif noseType.lower() == noseTypes[1]:  # Hemisphere
            diameter = find_parameter(file, "diameter")
            thickness = find_parameter(file, "thickness")
            density = find_parameter(file, "density")
            return Nose(noseType.lower(), eval(diameter), eval(thickness), eval(density))
        else:
            print("ERROR: invalid nose type '" + noseType + "' encountered at initialization. "
                                                            "Please check your spelling.")
            print("Possible options:")
            for string in noseTypes:
                print("\t" + string)
            exit(1)

class Body:
    def __init__(self, diameter, length, thickness, density):
        self.__diameter = diameter
        self.__length = length
        self.__density = density
        self.__thickness = thickness
        print("Body initialized!\n")

    def __str__(self):
        D = str(self.__diameter)
        l = str(self.__length)
        d = str(self.__thickness*1e3)
        m = str(self.getMass())
        rho = str(self.__density)
        return "Diameter: " + D + " m\n" + "Length: " + l + " m\n" + "Thickness: " + d + " mm\n" + "Mass: " + m + \
               " kg\n" + "Density: " + rho + " kgm^-3 "

    # Member functions
    def getSurfaceArea(self):
        R = self.__diameter/2
        L = self.__length
        return 2*np.pi*R*L

    def getVolume(self):
        d = self.__thickness
        R2 = self.__diameter/2
        R1 = R2 - d
        h = self.__length
        return np.pi*(R2**2 - R1**2)*h

    def getCavityVolum(self):
        d = self.__thickness
        R = self.__diameter/2 - d
        l = self.__length
        return np.pi*l*R**2

    def getLength(self):
        return self.__length

    def getDiameter(self):
        return self.__diameter

    def getMass(self):
        return self.__density*self.getVolume()

    def getInertiaMatrix(self):
        r2 = self.__diameter/2
        r1 = r2 - self.__thickness
        l = self.__length
        m = self.getMass()
        Ixx = 1/2*m*(r1**2 + r2**2)
        Iyy = 1/12*m*l**2
        Izz = Iyy
        return np.diag([Ixx, Iyy, Izz])

    def getCOM(self):
        l = self.__length
        return np.array([-l/2, 0, 0])  # COM relative to top of body

    @staticmethod
    def from_file(file):
        diameter = find_parameter(file, "diameter")
        length = find_parameter(file, "length")
        density = find_parameter(file, "density")
        thickness = find_parameter(file, "thickness")
        return Body(eval(diameter), eval(length), eval(thickness), eval(density))

class Fin:
    def __init__(self, *args):
        self.__semiChord = args[0]  # Fin semi span
        self.__rootChord = args[1]  # Fin root chord
        self.__tipChord = args[2]  # Fin tip chord
        self.__angle = args[3]  # Angle of ray from body to top outer edge
        self.__thickness = args[4]
        self.__density = args[5]
        print("Fin initialized!\n")

    def __str__(self):
        Chord = str(self.__semiChord)
        d = str(self.__thickness*1e3)
        m = str(self.getMass())
        rho = str(self.__density)
        return "Semi chord: " + Chord + " m\n" + "Thickness: " + d + " mm\n" + "Mass: " + m + " kg\n" \
               + "Density: " + rho + " kgm^-3"

    # Member functions
    def getSurfaceArea(self):
        l1 = self.__rootChord
        l2 = self.__tipChord
        cord = self.__semiChord
        return 1/2*(l1+l2)*cord

    def getVolume(self):
        l1 = self.__rootChord
        l2 = self.__tipChord
        cord = self.__semiChord
        d = self.__thickness
        return 1/2*(l1 + l2)*cord*d

    def getRootChord(self):
        return self.__rootChord

    def getTipChord(self):
        return self.__tipChord

    def getSemiChord(self):
        return self.__semiChord

    def getTopEdgeAngle(self):
        return self.__angle

    def getThickness(self):
        return self.__thickness

    def getMass(self):
        return self.__density*self.getVolume()

    def getInertiaMatrix(self):
        sigma = self.__density*self.__thickness  # Area density
        SC = self.__semiChord
        RC = self.__rootChord
        TC = self.__tipChord
        m = self.getMass()
        y = SC*(2*TC + RC)/(3*(TC + RC))
        Ixx = SC**3*(3*TC + RC)/12*sigma - m*y**2  # about root chord axis (along body edge)
        Iyy = Ixx/2
        Izz = Iyy
        return np.diag([Ixx, Iyy, Izz])

    def getCOM(self):
        cord = self.__semiChord
        l1 = self.__rootChord
        l2 = self.__tipChord
        a = self.__angle
        return np.array([-((cord/np.tan(a*np.pi/180) + (l1 + l2)/2)/2), 0, 0])  # COM relative to top edge of fin

    @staticmethod
    def from_file(file):
        semiChord = find_parameter(file, "semi_chord")
        rootChord = find_parameter(file, "root_chord")
        tipChord = find_parameter(file, "tip_chord")
        angle = find_parameter(file, "root_angle")
        density = find_parameter(file, "density")
        thickness = find_parameter(file, "thickness")
        return Fin(eval(semiChord), eval(rootChord), eval(tipChord), eval(angle), eval(thickness), eval(density))

class Motor:
    def __init__(self, *args):
        print("Initializing motor:")
        self.__name = args[0]
        self.__thrustMatrix = args[1]
        self.__timeArray = self.__thrustMatrix[:, 0]  # Assuming time values along 1st column
        self.__thrustArray = self.__thrustMatrix[:, 1]  # Assuming thrust values along 2nd column
        self.__diameter = args[3]
        self.__length = args[4]
        self.__initialPropellantMass = args[5]
        self.__frameMass = args[6]
        print("\tInterpolating thrust data...")
        self.__thrustFunction = interp1d(self.__timeArray, self.__thrustArray, kind='linear')  # Linear Interpolation for thrust curve
        self.__totalImpulse = args[2]
        self.__exhaustSpeed = self.__totalImpulse/self.__initialPropellantMass
        self.__burnTime = self.__timeArray[-1]
        self.__avgThrust = round(self.__totalImpulse/self.__burnTime, 4)
        dt = self.__burnTime/1e4  # This is usually in the order of 100 micro seconds (time step).
        timeList = np.arange(self.__timeArray[0], self.__burnTime + dt, dt)
        iterations = len(timeList)
        propellantMass = self.__initialPropellantMass
        massFlow = self.__thrustFunction(timeList)/self.__exhaustSpeed
        self.__propellantMassList = np.zeros(iterations)
        self.__propellantMassList[0] = propellantMass
        print("\tCalculating mass loss over the burn time of %1.2f s..." % self.__burnTime)
        for i in range(iterations - 1):
            propellantMass -= dt/2*(massFlow[i] + massFlow[i + 1])  # Trapezoid rule for integration of mdot over time
            self.__propellantMassList[i + 1] = propellantMass
        self.__propellantMassFunction = interp1d(timeList, self.__propellantMassList)
        print("Motor %s initialized!\n" % self.__name)

    def __str__(self):
        I = str(round(self.__totalImpulse, 2))
        avg = str(round(self.__avgThrust, 2))
        timeMax, Tmax = self.getMaxThrust()
        bTime = self.__burnTime
        name = self.getName()
        propMass = self.__initialPropellantMass*1e3
        frameMass = self.__frameMass*1e3
        sep = (len(name) + 10)*'-'
        return "Motor: " + name + "\n" + sep + "\nTotal impulse: " + I + " Ns\n" + "Average thrust: " + avg + " N" + \
               "\nMaximum thrust: " + str(Tmax) + " N,\tat time " + str(timeMax) + " s\n" + "Burntime: " + \
               str(bTime) + " s\n" + "Propellant mass: " + str(propMass) + " g\n" + "Frame mass: " + str(frameMass) + " g\n"

    # Get functions
    def getName(self):
        return self.__name

    def getPropellantMass(self, t):
        if t <= self.__timeArray[0]:
            return self.__initialPropellantMass
        elif t >= self.__timeArray[-1]:
            return self.__propellantMassList[-1]

        return self.__propellantMassFunction(t)

    def getMass(self, t):
        return self.getPropellantMass(t) + self.__frameMass

    def getInertiaMatrix(self, t):
        """
        return the moment of inertia about principle axes with origin at COM
        """
        frameMass = self.__frameMass
        if t >= self.__timeArray[-1]:
            propMass = 0
        elif t < self.__timeArray[0]:
            propMass = self.__initialPropellantMass
        else:
            propMass = self.__propellantMassFunction(t)

        r0 = self.getDiameter()/2
        h = self.getLength()
        fuelMassRatio = propMass/self.__initialPropellantMass
        Ixx = frameMass*r0**2 + propMass*(r0**2)/2*(2 - fuelMassRatio)
        Iyy = propMass*((h**2)/12 + 1/4*(r0**2)*(2 - fuelMassRatio)) + frameMass*((h**2)/12 + 1/4*(r0**2))
        Izz = Iyy
        return np.diag([Ixx, Iyy, Izz])
        #TODO implement MOI for motor with axial fuel burn.

    def getCOM(self, t, fuelBurnRadially=True):
        l = self.__length
        # If fuel burns radially, quickly return the well know solution (relative to top of motor)
        if fuelBurnRadially:
            return np.array([-l/2, 0, 0])

        frameMass = self.__frameMass
        propMass = self.__initialPropellantMass
        if t < self.__timeArray[0]:
            Mtot = propMass + frameMass
        elif t >= self.__timeArray[-1]:
            propMass = self.__propellantMassList[-1]
            Mtot = frameMass + propMass
        else:
            propMass = self.__propellantMassFunction(t)
            Mtot = propMass + frameMass

        COM = -1/Mtot*(frameMass + propMass**2/self.__initialPropellantMass)*l/2  # COM relative to top of motor
        return np.array([COM, 0, 0])

    def getLength(self):
        return self.__length

    def getDiameter(self):
        return self.__diameter

    def getMaxThrust(self):
        i = np.argmax(self.__thrustArray)
        maxThrust = self.__thrustMatrix[i]
        return maxThrust[0], maxThrust[1]  # 1 timeAtMax, 2 maxThrust

    def getAverageThrust(self):
        return self.__avgThrust

    def getTotalImpulse(self):
        return self.__totalImpulse

    def getBurnTime(self):
        return self.__burnTime

    # Auxiliary functions
    def thrust(self, t):
        if t <= self.__burnTime:
            return self.__thrustFunction(t)
        else:
            return 0

    def massFlow(self, t):
        return self.thrust(t)/self.__exhaustSpeed

    def plotPerformance(self, show=True):
        dt = self.__burnTime/1e4
        timeList = np.arange(self.__timeArray[0], self.__burnTime + dt, dt)
        thrustArray = self.__thrustFunction(timeList)
        propellantMassArray = self.__propellantMassList
        COMarray = np.array([self.getCOM(t)[0] for t in timeList])
        # PLOT FORCE
        plt.figure()
        ax1 = plt.subplot(211, xlabel='time [s]', ylabel='[N]')
        ax1.plot(timeList, thrustArray, label='Thrust', c='r', lw='2')
        ax1.set_title(r'Thrust during burn phase of %s' %self.__name)
        ax1.grid()
        ax1.legend(loc='best')
        # PLOT PROPELLANT MASS LOSS
        ax2 = plt.subplot(212, xlabel='time [s]', ylabel='[g]')
        ax2.plot(timeList, propellantMassArray*1e3, label='propellant mass', c='b', lw='2')
        ax2.set_title('Propellant mass')
        ax2.grid()
        ax2.legend(loc='best')
        plt.subplots_adjust(hspace=0.5)
        # PLOT COM OVER TIME
        plt.figure()
        plt.plot(timeList, COMarray*100, label='COM', c='r', lw='2')
        plt.title('COM of %s during burn phase, length %1.1f cm' % (self.__name, self.__length*100))
        plt.grid()
        plt.legend(loc='best')
        if show:
            plt.show()

    @staticmethod
    def from_file(motorFile):
        """
            Read a file with motor specs.
            ASSUMPTIONS: -
            :param motorFile: The file with motor specs in proper format
            :return: motor
            """
        name = find_parameter(motorFile, "name")
        diameter = find_parameter(motorFile, "diameter")
        length = find_parameter(motorFile, "length")
        propMass = find_parameter(motorFile, "propellant_mass")
        frameMass = find_parameter(motorFile, "frame_mass")
        totalImpulse = find_parameter(motorFile, "total_impulse")
        motorFile = open(motorFile, 'r')

        thrust = [[0, 0]]
        with motorFile as fp:
            for i, line in enumerate(fp):
                if i > 5:
                    row = line.split()
                    t = eval(row[0])
                    f = eval(row[1])
                    thrust.append([t, f])
        thrust = np.array(thrust)
        return Motor(name, thrust, float(totalImpulse), float(diameter), float(length), float(propMass), float(frameMass))

class Payload:
    def __init__(self, width):
        self.__mass = 4  # this mass is fixed for all rockets qualified for competition.
        self.__width = width
        print('Payload initialized!\n')

    def getMass(self):
        return self.__mass

    def getInertiaMatrix(self):
        m = self.__mass
        r = self.__width/2
        Ixx = 2/5*m*r**2  # Assuming rough spherical shape
        Iyy = Ixx
        Izz = Iyy
        return np.diag([Ixx, Iyy, Izz])

    def getWidth(self):
        return self.__width

    @staticmethod
    def from_file(file=''):
        width = find_parameter(file, "width")
        return Payload(eval(width))

class RocketSimple:
    def __init__(self, nose, body, fin, numberOfFins, motor, payload, partsPlacement):
        print("Initalizing rocket:")
        self.__partsPlacement = partsPlacement
        self.__rocketStructure = np.array([nose, payload, body, fin])
        self.__rocketMotor = motor
        self.__N = numberOfFins  # Number of fins on rocket
        self.__massOfRocketStructure = np.array([part.getMass() for part in self.__rocketStructure])
        self.__massOfRocketStructure[3] = self.__N*self.__massOfRocketStructure[3]  # There are N fins
        print("\tCalculating rocket mass..")
        # (add 4 kg for now to account for electronics/recovery etc.)
        # TODO Account for electronics/recovery etc.
        self.__rocketMass = self.__massOfRocketStructure.sum() + 4
        self.__motorMass = self.__rocketMotor.getMass(0)
        # TOTAL MASS
        self.__mass = self.__rocketMass + self.__motorMass

        # COM
        print("\tCalculating rocket COM (relative to rocket origin)..")
        self.__noseCOM = nose.getCOM()[0] - nose.getLength()
        self.__bodyCOM = body.getCOM()[0] - nose.getLength()
        # Assuming placement of fin is position of bottom edge relative to body top
        self.__finCOM = partsPlacement[0] + (fin.getRootChord() + fin.getCOM()[0]) - nose.getLength()
        # Assuming placement of motor is at bottom of rocket (motor bottom align with body bottom)
        self.__motorCOM = motor.getLength() + motor.getCOM(0)[0] - body.getLength() - nose.getLength()
        # Assuming placement of payload is its COM relative to body top
        self.__payloadCOM = partsPlacement[1] - nose.getLength()
        self.__COMofRocketStructure = np.array([self.__noseCOM, self.__payloadCOM, self.__bodyCOM, self.__finCOM])
        self.__rocketStructureCOM = np.dot(self.__massOfRocketStructure, self.__COMofRocketStructure)/self.__rocketMass
        # FINAL COM OF ROCKET
        self.__COM = (self.__rocketStructureCOM*self.__rocketMass + self.__motorCOM*self.__motorMass)/self.__mass

        # Nose COP
        self.__Xnose = 0
        self.__CNnose = 2
        if nose.getNoseType() == noseTypes[0]:  # Conic
            self.__Xnose = -0.666*nose.getLength()
        elif nose.getNoseType() == noseTypes[1] or nose.getNoseType() == noseTypes[2]:  # Hemisphere or Ogive
            self.__Xnose = -0.446*nose.getLength()
        self.__Xcp_nose = self.__Xnose

        # Body COP
        Aref = np.pi*(body.getDiameter()/2)**2  # Crossection of body
        Aplan_body = body.getDiameter()*body.getLength()
        Aplan_nose = 2/3*nose.getLength()*nose.getDiameter()
        Xplan_body = nose.getLength() + 1/2*body.getLength()
        Xplan_nose = 5/8*nose.getLength()
        Aplan_total = Aplan_body + Aplan_nose
        self.__Xcp_body = -(Xplan_body*Aplan_body + Xplan_nose*Aplan_nose)/Aplan_total
        self.__CNbody = Aplan_total/Aref

        # Fin COP
        R = body.getDiameter()/2  # Radius of body
        SC = fin.getSemiChord()  # Semi chord of fins
        RC = fin.getRootChord()
        TC = fin.getTipChord()
        theta = fin.getTopEdgeAngle()*np.pi/180.0 # TopEdgeAngle is in degrees
        y = SC*(2*TC + RC)/(3*(TC + RC))
        Lf = np.sqrt(SC**2 + (SC/np.tan(theta) + 1/2*(TC - RC))**2)
        self.__CNfin = (1 + R/(R + SC))*(4*self.__N*(SC/(2*R))**2/(1 + np.sqrt(1 + (2*Lf/(RC + TC))**2)))
        Xb = partsPlacement[0] - nose.getLength()
        Xr = -SC/np.tan(theta)
        Xf = Xb + Xr/3*(RC + 2*TC)/(RC + TC) - 1/6*((RC + TC) - RC*TC/(RC + TC))
        self.__Xcp_fin = Xf

        print("\tCalculating inertia matrix of rocket..")
        # MOMENT OF INERTIA (about rocket axes with origin at COM, calculated with parallel axis thm)
        noseMOI = nose.getInertiaMatrix() + np.diag([0, 1, 1])*nose.getMass()*(self.__COM - self.__noseCOM)**2
        bodyMOI = body.getInertiaMatrix() + np.diag([0, 1, 1])*body.getMass()*(self.__COM - self.__bodyCOM)**2
        finMOI = self.__N*(fin.getInertiaMatrix() + (np.diag([1, 0, 0])*(body.getDiameter()/2 + y)**2 +
                                                     np.diag([0, 1, 1])*(self.__COM - self.__finCOM)**2)*fin.getMass())
        motorMOI = motor.getInertiaMatrix(0) + np.diag([0, 1, 1])*self.__motorMass*(self.__COM - self.__motorCOM)**2
        payloadMOI = payload.getInertiaMatrix() + np.diag([0, 1, 1])*payload.getMass()*(self.__COM - self.__payloadCOM)**2

        self.__rocketStructureMOI = noseMOI + bodyMOI + finMOI + payloadMOI
        # FINAL MOMENT OF INERTIA MATRIX
        self.__InertiaMatrix = self.__rocketStructureMOI + motorMOI

        # TOTAL LENGTH OF ROCKET
        self.__length = nose.getLength() + body.getLength() + (SC/np.tan(theta)-RC) + TC
        # MAXIMAL WIDTH OF ROCKET
        self.__width = body.getDiameter()/2 + SC
        # DRAG COEFFICIENT (at some arbitrary speed, 30 m/s)
        Forces.updateCd_2(self, [0, 0, 0], [30, 0, 0], 0)
        print("Rocket initialized!\n")
        self.printSpecifications(0, 5*np.pi/180) # Specs at AoA = 5 deg.

    # Rocket parts
    def getNose(self):
        return self.__rocketStructure[0]

    def getPayload(self):
        return self.__rocketStructure[1]

    def getBody(self):
        return self.__rocketStructure[2]

    def getFin(self):
        return self.__rocketStructure[3]

    def getMotor(self):
        return self.__rocketMotor

    # Rocket structure
    def getNumberOfFins(self):
        return self.__N

    def getMass(self, t):
        self.__mass = self.__rocketMass + self.__rocketMotor.getMass(t)
        return self.__mass

    def setMass(self, mass):
        self.__rocketMass = mass

    def getInertiaMatrix(self, t):
        self.__motorCOM = self.__rocketMotor.getLength() + self.__rocketMotor.getCOM(t)[0] - self.__rocketStructure[
            0].getLength() - self.__rocketStructure[2].getLength()
        motorMass = self.getMotor().getMass(t)
        motorMOI = self.getMotor().getInertiaMatrix(t) + np.diag([0, 1, 1])*motorMass*(self.getCOM(t)[0] - self.__motorCOM)**2
        self.__InertiaMatrix = self.__rocketStructureMOI + motorMOI
        return self.__InertiaMatrix

    def getLength(self):
        return self.__length

    def getTransversalExtension(self):
        "The distance from center of body to fin tip"
        return self.__width

    def getCOM(self, t):
        mass = self.getMass(t)
        self.__motorCOM = self.__rocketMotor.getLength() + self.__rocketMotor.getCOM(t)[0] - self.__rocketStructure[
            0].getLength() - self.__rocketStructure[2].getLength()
        motorMass = self.getMotor().getMass(t)
        self.__COM = (self.__rocketStructureCOM*self.__rocketMass + self.__motorCOM*motorMass)/mass
        return np.array([self.__COM, 0, 0])

    def setCOM(self, com):
        self.__COM = com

    def getCOMofParts(self):
        return self.__COMofRocketStructure

    # Aerodynamics
    def getCOP(self, AoA):
        """
        :param AoA: [float] the angle of attack [rad]

        :return: [np.array] Position of COP relative to nose tip
        """
        #Nose, Body and Fins Cn
        CNnose = self.__CNnose*np.sin(AoA)
        CNfin = self.__CNfin*AoA
        CNbody = self.__CNbody*np.sin(AoA)**2
        CNrocket = CNnose + CNbody + CNfin
        COP0 = (CNnose*self.__Xcp_nose + CNbody*self.__Xcp_body + CNfin*self.__Xcp_fin)/CNrocket
        return np.array([COP0, 0, 0])

    def getStabilityMargin(self, AoA, t=0):
        COM = self.getCOM(t)[0]
        COP = self.getCOP(AoA)[0]
        return COM - COP

    def getAeroForces(self, AoA, position, velocity):
        """
        :param velocity: [np.array] velocity of rocket (with wind) relative to world [m/s]
        :param AoA: [float] the angle of attack [rad]

        :return: [np.array] ([drag, lift]) on rocket attacking in COP [N]
        """
        Forces.updateCd_2(self, position, velocity, AoA)
        drag = Forces.SAMdrag(self, position, velocity)
        lift = Forces.SAMlift(self, position, velocity, AoA)
        return np.array([drag, lift])

    def getMomentAboutCOM(self, position, velocity, AoA):
        """
        :param velocity: [float] the air speed relative to rocket [m/s]
        :param AoA: [float] the angle of attack [rad]

        :return: [np.array] The total moment on rocket about COM (component normal to aerodynamic plane) [Nm]
        """
        #TODO Implement this

    def getCd(self):
        return self.__Cd

    def getCn(self, AoA):
        #Nose, Body and Fins Cn
        CNnose = self.__CNnose*np.sin(AoA)
        CNfin = self.__CNfin*AoA
        CNbody = self.__CNbody*np.sin(AoA)**2
        # FINAL ROCKET COP
        CNrocket = CNnose + CNbody + CNfin
        return CNrocket

    # Set functions
    def setCd(self, Cd):
        self.__Cd = Cd

    # auxiliary
    def printSpecifications(self, t, AoA=0):
        Mass = self.getMass(t)
        motorCOM = self.__motorCOM + self.getMotor().getCOM(t) - self.getMotor().getCOM(0)
        COM = self.getCOM(t)[0]
        COP = self.getCOP(AoA)[0]
        stability_margin = self.getStabilityMargin(AoA)/self.getBody().getDiameter()
        MOI = self.getInertiaMatrix(t)
        Cd = self.getCd()
        Cn = self.getCn(AoA)
        partsCOM = np.append(self.getCOMofParts(), motorCOM)
        partsMass = self.__massOfRocketStructure
        partsMass = np.append(partsMass, self.__motorMass)
        partNames = ['Nose', 'Payload','Body', 'Fins', 'Motor']
        length = self.getLength()
        max_extension = self.getTransversalExtension()
        N = self.__N  # Number of fins
        dots = 33*'-'

        print(dots)
        print("Rocket Specifications at time %1.1f" % t)
        print(dots)
        print("Mass: %1.2f kg" % Mass)
        print("Moment of inertia (about rocket axes with COM as origin) [kg*cm^2]:")
        print(np.array2string(MOI*1e4, precision=1))
        print("Length: %1.1f cm" % (length*1e2))
        print("Transversal extension: %1.2f cm" % (max_extension*1e2))
        print("Number of fins: %d" % N)
        print("Drag coeff: %1.2f" % Cd)
        print("Lift coeff: %1.2f (at AoA=%1.1f deg)" % (Cn, AoA*180/np.pi))
        print("COM of (x-coordinate): ")
        for i in range(len(partNames)):
            print("\t%s: %1.1f cm (%1.2f kg)" % (partNames[i], partsCOM[i]*1e2, partsMass[i]))
        print("COM (x-coordinate): %1.1f cm" % (COM*1e2))
        print("COP (x-coordinate): %1.1f cm" % (COP*1e2))
        print("Stability margin: %1.2f b.c" % (stability_margin)) # In body calibers
        print(dots)

    @staticmethod
    def from_file(rocket_file, path_to_file=""):
        """
        Creating an instance of a rocket by reading a rocket file that is located in a folder containing files for all
        necessary rocket parts.
        Example: folder 'myRocket'
                    'myRocket' content:
                                - myNose.dot
                                - myBody.dot
                                - myMotor.dot
                                - myFin.dot
                                - myPayload.dot
                                - myRocket.dot   <--  rocket file!
                 If this folder is located in the current working folder, you can simply create an instance of this
                 rocket with:
                          myRocket = RocketSimple.from_file('myRocket.dot', 'myRocket')
        :param path_to_file: [string] a path to the rocket file relative to the current path (empty by default)
        :param rocket_file: [string] name of rocket file
        :return: [RocketSimple class] Rocket instance with specs from rocket file.
        """
        # get file names of each rocket part
        path = path_to_file + rocket_file
        noseFile = find_parameter(path, "nose")
        bodyFile = find_parameter(path, "body")
        finFile = find_parameter(path, "fin")
        numberOfFins = find_parameter(path, "number_of_fins")
        motorFile = find_parameter(path, "motor")
        payloadFile = find_parameter(path, "payload")
        payloadPlacement = find_parameter(path, "payload_placement")
        finPlacement = find_parameter(path, "fin_placement")
        partsPlacement = np.array([eval(finPlacement), eval(payloadPlacement)])

        # Initialize rocket parts
        path = path_to_file
        nose = Nose.from_file(path + noseFile)
        body = Body.from_file(path + bodyFile)
        fin = Fin.from_file(path + finFile)
        motor = Motor.from_file(path + motorFile)
        payload = Payload.from_file(path + payloadFile)

        return RocketSimple(nose, body, fin, eval(numberOfFins), motor, payload, partsPlacement)
