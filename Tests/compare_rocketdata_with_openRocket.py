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

sample_file = 'full_report_edited.dot'
init_file = 'mass_properties_rocket_v9_edited.dot'
path = 'V9/'
rocket = Rocket.from_file_with_AoAspeed(init_file, sample_file, path)
BT = rocket.getMotor().getBurnTime()

openRocket_file = 'V9/openRocketV9data.txt'
OR_time, OR_mass, OR_propMass, OR_COM, OR_Ixx, OR_Ilong = unwrap_openRocket(openRocket_file)

n = len(OR_time)
r_mass = np.zeros(n)
r_propMass = np.zeros(n)
r_COM = np.zeros(n)
r_Ixx = np.zeros(n)
r_Ilong = np.zeros(n)
for i in range(n):
        r_mass[i] = rocket.getMass(OR_time[i])*1e3  # In grams
        r_propMass[i] = rocket.getMotor().getPropellantMass(OR_time[i])*1e3  # In grams
        r_COM[i] = rocket.getCOM(OR_time[i])*1e2  # In cm
        I = rocket.getInertiaMatrix(OR_time[i])  # In kgm^2
        r_Ixx[i] = I[0][0]  # In kgm^2
        r_Ilong[i] = I[1][1]  # In kgm^2

print(r_mass)
print(OR_mass)

# PLotting
plt.figure()
ax1 = plt.subplot(211, xlabel='time [s]', ylabel='[g]')
ax1.plot(OR_time, r_mass - OR_mass, label='Delta_Mass(t)', lw=2, c='r')
ax1.set_title('Mass difference Propulse Rocket and OpenR. during burntime of %1.1f s'%BT)
ax1.grid()

ax2 = plt.subplot(212, xlabel='time [s]', ylabel='[g]', sharex=ax1)
ax2.plot(OR_time, r_propMass - OR_propMass, label='Delta_propellantMass(t)', lw=2, c='b')
ax2.set_title('Propellant mass difference')
ax2.grid()
plt.subplots_adjust(hspace=0.5)

plt.figure()
ax3 = plt.subplot(211, xlabel='time [s]', ylabel='[cm]')
ax3.plot(OR_time, -r_COM - OR_COM, label='Delta_COM(t)', lw=2, c='r')
ax3.set_title('COM of rocket during burn phase of %1.1f s'%BT)
ax3.grid()

plt.figure()
ax4.figure()




