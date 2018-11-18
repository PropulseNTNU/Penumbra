import sys
sys.path.append('../Rocket/')
sys.path.append('../Forces/')
sys.path.append('../Trajectory/')
import numpy as np
import Forces as Forces
import Trajectory as Trajectory
from Rocket2 import Rocket
import matplotlib.pyplot as plt

def test_trajectory_module():
    rad2deg = 180/np.pi
    sample_file = 'full_report_edited.dot'
    init_file = 'mass_properties_rocket_v9_edited.dot'
    path = 'V9/'
    rocket = Rocket.from_file_with_AoAspeed(init_file, sample_file, path)
    initialInclination = 6.0/180.0*np.pi
    launchRampLength = 2*rocket.getLength()
    timeStep = 0.005
    simulationTime= 7
    (t, position, euler, linearVelocity, angularVelocity, AoA, thrust, gravity, drag, lift) \
    = Trajectory.calculateTrajectory(rocket, initialInclination, launchRampLength, timeStep, simulationTime)

    plt.figure()
    ax1 = plt.subplot(311, xlabel='time [s]', ylabel='x [m]')
    ax1.plot(t, position[:,0], label='x', lw=2, c='r')
    ax1.set_title('x')
    ax1.grid()
    plt.subplots_adjust(hspace=0.5)
    ax2 = plt.subplot(312, xlabel='time [s]', ylabel='y [m]')
    ax2.plot(t, position[:,1], label='y', lw=2, c='b')
    ax2.set_title('y')
    ax2.grid()
    plt.subplots_adjust(hspace=0.5)
    ax3 = plt.subplot(313, xlabel='time [s]', ylabel=' h [m]')
    ax3.plot(t, -position[:,2], label='h', lw=2, c='g')
    ax3.set_title('h')
    ax3.grid()

    plt.figure()
    ax1 = plt.subplot(111, xlabel='r [m]', ylabel=' h [m]')
    ax1.plot(np.sqrt(position[:,0]**2+position[:,1]**2), -position[:,2], lw=2, c='r')
    ax1.set_title('h vs. r')
    ax1.grid()
    ax1.axis('equal')

    plt.figure()
    ax1 = plt.subplot(411, xlabel='time [s]', ylabel='pitch [deg]')
    ax1.plot(t, euler[:,0]*rad2deg, label='pitch', lw=2, c='r')
    ax1.set_title('pitch')
    ax1.grid()
    plt.subplots_adjust(hspace=0.5)
    ax2 = plt.subplot(412, xlabel='time [s]', ylabel='yaw [deg]')
    ax2.plot(t, euler[:,1]*rad2deg, label='yaw', lw=2, c='b')
    ax2.set_title('yaw')
    ax2.grid()
    plt.subplots_adjust(hspace=0.5)
    ax3 = plt.subplot(413, xlabel='time [s]', ylabel='roll [deg]')
    ax3.plot(t, euler[:,2]*rad2deg, label='roll', lw=2, c='g')
    ax3.set_title('roll')
    ax3.grid()
    plt.subplots_adjust(hspace=0.5)
    ax1 = plt.subplot(414, xlabel='time [s]', ylabel='AoA [deg]')
    ax1.plot(t[0:-1], AoA*rad2deg, label='AoA', lw=2, c='k')
    ax1.set_title('AoA')
    ax1.grid()

    plt.figure()
    ax1 = plt.subplot(311, xlabel='time [s]', ylabel='v_x [m/s]')
    ax1.plot(t, linearVelocity[:,0], label='COMx(t)', lw=2, c='r')
    ax1.set_title('v_x')
    ax1.grid()
    plt.subplots_adjust(hspace=0.5)
    ax2 = plt.subplot(312, xlabel='time [s]', ylabel='v_y [m/s]')
    ax2.plot(t, linearVelocity[:,1], label='mass(t)', lw=2, c='b')
    ax2.set_title('v_y')
    ax2.grid()
    plt.subplots_adjust(hspace=0.5)
    ax3 = plt.subplot(313, xlabel='time [s]', ylabel='v_h [m/s]')
    ax3.plot(t, -linearVelocity[:,2], label='mass(t)', lw=2, c='g')
    ax3.set_title('v_h')
    ax3.grid()

    plt.figure()
    ax1 = plt.subplot(311, xlabel='time [s]', ylabel='w_x [rad/s]')
    ax1.plot(t, angularVelocity[:,0], label='COMx(t)', lw=2, c='r')
    ax1.set_title('w_x')
    ax1.grid()
    plt.subplots_adjust(hspace=0.5)
    ax2 = plt.subplot(312, xlabel='time [s]', ylabel='w_y [rad/s]')
    ax2.plot(t, angularVelocity[:,1], label='mass(t)', lw=2, c='b')
    ax2.set_title('w_y')
    ax2.grid()
    plt.subplots_adjust(hspace=0.5)
    ax3 = plt.subplot(313, xlabel='time [s]', ylabel='w_z [rad/s]')
    ax3.plot(t, angularVelocity[:,2], label='mass(t)', lw=2, c='g')
    ax3.set_title('w_z')
    ax3.grid()

    plt.figure()
    ax1 = plt.subplot(311, xlabel='time [s]', ylabel='[N]')
    ax1.plot(t[0:-1], thrust[:,0], label='thrust', lw=2, c='r')
    ax1.plot(t[0:-1], drag[:,0], label='drag', lw=2, c='b')
    ax1.plot(t[0:-1], lift[:,0], label='lift', lw=2, c='g')
    ax1.plot(t[0:-1], gravity[:,0], label='gravity', lw=2, c='k')
    ax1.set_title('Forces x')
    ax1.grid()
    plt.subplots_adjust(hspace=0.5)
    ax2 = plt.subplot(312, xlabel='time [s]', ylabel='[N]')
    ax2.plot(t[0:-1], thrust[:,0], label='thrust', lw=2, c='r')
    ax2.plot(t[0:-1], drag[:,1], label='drag', lw=2, c='b')
    ax2.plot(t[0:-1], lift[:,1], label='lift', lw=2, c='g')
    ax2.plot(t[0:-1], gravity[:,1], label='gravity', lw=2, c='k')
    ax2.set_title('Forces y')
    ax2.grid()
    plt.subplots_adjust(hspace=0.5)
    ax3 = plt.subplot(313, xlabel='time [s]', ylabel='[N]')
    ax3.plot(t[0:-1], thrust[:,2], label='thrust', lw=2, c='r')
    ax3.plot(t[0:-1], drag[:,2], label='drag', lw=2, c='b')
    ax3.plot(t[0:-1], lift[:,2], label='lift', lw=2, c='g')
    ax3.plot(t[0:-1], gravity[:,2], label='gravity', lw=2, c='k')
    ax3.set_title('Forces z')
    ax3.grid()

    plt.show()

def main():
    test_trajectory_module()

main()
