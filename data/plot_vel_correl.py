import matplotlib as matl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pylab as P
from math import *

vel_correl = {}

# Data for gas
data = np.loadtxt("vel_corrGas.dat", skiprows=5)
vel_correl['GAS'] = {}
vel_correl['GAS']['VEC'] = data[:,0]
vel_correl['GAS']['MOD'] = data[:,2]

# Data for liquid
data = np.loadtxt("vel_corrLiquid.dat", skiprows=5)
vel_correl['LIQ'] = {}
vel_correl['LIQ']['VEC'] = data[:,0]
vel_correl['LIQ']['MOD'] = data[:,2]

# Data for solid
data = np.loadtxt("vel_corrSolid.dat", skiprows=5)
vel_correl['SOL'] = {}
vel_correl['SOL']['VEC'] = data[:,0]
vel_correl['SOL']['MOD'] = data[:,2]

# Vector correlation
f, axarr = plt.subplots(3, sharex=True)
indexes = range(0,400)
x = [0.004*y for y in indexes]

axarr[0].plot(x, vel_correl['GAS']['VEC'][indexes], 'k-')
axarr[0].text(0.3, 0.6, 'Gas', color = 'black')
axarr[1].plot(x, vel_correl['LIQ']['VEC'][indexes], 'k-')
axarr[1].text(0.3, 0.6, 'Liquid', color = 'black')
axarr[2].plot(x, vel_correl['SOL']['VEC'][indexes], 'k-')
axarr[2].text(0.3, 0.6, 'Solid', color = 'black')
axarr[2].set_xlabel(r'Time ($\tau$)', fontsize=15)

# Modulus correlation
f, axarr = plt.subplots(2, sharex=True)
indexes = range(0,500)
x = [0.004*y for y in indexes]

axarr[0].loglog(x, vel_correl['GAS']['MOD'][indexes], 'k-')
axarr[0].text(0.3, 0.5, 'Gas', color = 'black')

indexes = range(0,200)
x = [0.004*y for y in indexes]
index_lin = range(30,150)
t_lin = [0.004*y for y in index_lin]
z = np.polyfit(np.log(t_lin),np.log(vel_correl['LIQ']['MOD'][index_lin]),1)
t_plot = [0.004*y for y in range(30,200)]
fitted = np.exp([z[1] - 1.5*log(y) for y in t_plot])
print z
axarr[1].loglog(x, vel_correl['LIQ']['MOD'][indexes], 'k--')
axarr[1].text(0.1, 0.5, 'Liquid', color = 'black')
axarr[1].loglog(t_plot, fitted, 'r-')
axarr[1].text(0.9, 0.01, '3/2 decay', color = 'red')
axarr[1].set_xlabel(r'Time ($\tau$)', fontsize=15)

plt.show()