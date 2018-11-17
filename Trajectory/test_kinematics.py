import numpy as np
import Kinematics

def test():
    deg2rad = np.pi/180
    rad2deg = 180/np.pi
    R0 = Kinematics.Ryzx(0,0,0)
    R1 = Kinematics.Ryzx(np.pi/4,0,0)
    R2 = Kinematics.Ryzx(0,np.pi/4,0)
    R3 = Kinematics.Ryzx(0,0,np.pi/4)
    print(R0)
    print(R1)
    print(R2)
    print(R3)
    R4 = Kinematics.Ryzx(-80*deg2rad,0,0)
    print(R4)
    pitch = np.pi/3
    yaw = -np.pi/7
    roll = np.pi/4
    R5 = Kinematics.Ryzx(pitch, yaw, roll)
    q = Kinematics.euler2quaternion(pitch, yaw, roll)
    R6 = Kinematics.Rquaternion(q)
    print(R5)
    print(R6)
    print(R5-R6)
    w = np.array([1.5, -2.3, 3.1])
    print(Kinematics.CrossProductMatrix(w))
    print(np.dot(Kinematics.CrossProductMatrix(w),w.T))
    H = Kinematics.TransformationMatrix(w.T)
    print(H)

def main():
    test()

main()
