#!/usr/bin/python

# Declare libraries
import matplotlib as matl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pylab as P
from math import *

LINES_HEADER = 5
gas = np.loadtxt("correlationGas.dat", skiprows=LINES_HEADER)
solid = np.loadtxt("correlationSolid.dat", skiprows=LINES_HEADER)
liquid = np.loadtxt("correlationLiquid.dat", skiprows=LINES_HEADER)

f, axarr = plt.subplots(2, sharex=True)
indexes = range(0,400)
x = [0.01*y for y in indexes]
axarr[0].plot(x, gas[indexes], 'r-')
axarr[0].plot(x, liquid[indexes], 'k--')
axarr[0].text(1.2, 2.5, 'Liquid', color = 'black')
axarr[0].text(1.5, 1.25, 'Gas', color = 'red')
axarr[1].plot(x, solid[indexes], 'k-')
axarr[1].text(1.1, 7, 'Solid', color = 'black')
axarr[1].set_xlabel(r'$r/\sigma$', fontsize=15)

plt.show()