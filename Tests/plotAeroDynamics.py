"""
Plots the interpolation of Forces/moment

Last edit: 17.11.18
"""
import sys
sys.path.append('../Rocket/')
from Rocket2 import Rocket

from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import numpy as np


sample_file = 'full_report_edited.dot'
init_file = 'mass_properties_rocket_v9_edited.dot'
path = 'V9/'
rocket = Rocket.from_file_with_AoAspeed(init_file, sample_file, path)

speed = np.arange(0, 0.5, 5e-3)
aoa = np.arange(0, 10, 5e-3)

# Aero data
angle, vel, drag, lift, moment = rocket.getAeroData()

# Forces and moment
Lift = rocket.getLift(aoa, speed)
Drag = rocket.getDrag(aoa, speed)
Moment = rocket.getMomentAboutCOM(aoa, speed)
AoA, Speed = np.meshgrid(aoa, speed)

# Surface
fig = plt.figure()
ax = fig.gca(projection='3d')

ax.scatter3D(vel, angle, lift)
ax.plot_surface(Speed, AoA, Lift, cmap=cm.coolwarm,
                      lw=0, antialiased=False)

ax.set_xlabel('Speed [mach]')
ax.set_ylabel('AoA [degree]')
ax.set_zlabel('Lift [N]')
plt.show()
