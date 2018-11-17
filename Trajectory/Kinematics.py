import numpy as np

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
