"""

This program compares the time evolution of the rocket in propulse simulator with
the OpenRocket simulator.

Last edit: 18.11.2018

"""
import sys
sys.path.append('../Rocket/')
from Rocket2 import Rocket
from lib.File_utilities import unwrap_openRocket

import numpy as np
import matplotlib.pyplot as plt

sample_file = 'V9_CFD.txt'
init_file = 'V9_data.dot'
path = 'V9/'
rocket = Rocket.from_file_with_AoAspeed(init_file, sample_file, path)
BT = rocket.getMotor().getBurnTime()

openRocket_file = 'V9/openRocketV9data.txt'
OR_time, OR_mass, OR_propMass, OR_COM, OR_Ixx, OR_Ilong = unwrap_openRocket(openRocket_file)

n = len(OR_time)
r_mass = np.zeros(n)
r_COM = np.zeros(n)
r_Ixx = np.zeros(n)
r_Ilong = np.zeros(n)
for i in range(n):
        r_mass[i] = rocket.getMass(OR_time[i])*1e3  # In grams
        r_COM[i] = rocket.getCOMx(OR_time[i])*1e2  # In cm
        I = rocket.getInertiaMatrix(OR_time[i])  # In kgm^2
        r_Ixx[i] = I[0][0]*1e4  # In kgm^2
        r_Ilong[i] = I[1][1]*1e4  # In kgm^2

# PLotting
plt.figure()
ax1 = plt.subplot(111, xlabel='time [s]', ylabel='[g]')
ax1.plot(OR_time, (r_mass - OR_mass), label='dMass(t)', lw=2, c='r')
ax1.set_title('Absolute error in Propulse Rocket and OpenR. mass during burntime of %1.1f s'%BT)
ax1.legend(loc='best')
ax1.grid()

plt.figure()
plt.plot(OR_time, (r_COM + OR_COM)/r_COM*100, label='dCOM(t)', lw=2, c='g')
plt.xlabel('time [s]')
plt.ylabel('[%]')
plt.title('Relative error of COM during burntime of %1.1f s'%BT)
plt.grid()
plt.legend(loc='best')

plt.figure()
ax4 = plt.subplot(211, xlabel='time [s]', ylabel='[%]')
ax4.plot(OR_time, (r_Ixx - OR_Ixx*1e4)/r_Ixx*100, label='dI_roll (t)', lw=2, c='b')
ax4.set_title('Relative error for MOI during burntime of %1.1f s'%BT)
ax4.legend(loc='best')
ax4.grid()

ax5 = plt.subplot(212, xlabel='time [s]', ylabel='[%]', sharex=ax4)
ax5.plot(OR_time, (r_Ilong - OR_Ilong*1e4)/r_Ilong*100, label='dI_longi. (t)', lw=2, c='b')
ax5.legend(loc='best')
ax5.grid()

# Plot settings
plt.subplots_adjust(hspace=0.5)
plt.show()
