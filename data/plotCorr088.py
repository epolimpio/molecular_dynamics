#!/usr/bin/python

# Declare libraries
import matplotlib as matl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pylab as P
from math import *

LINES_HEADER = 5
cor1 = np.loadtxt("correlation013.dat", skiprows=LINES_HEADER)
cor2 = np.loadtxt("correlation021.dat", skiprows=LINES_HEADER)
cor3 = np.loadtxt("correlation005.dat", skiprows=LINES_HEADER)

indexes = range(0,400)
x = [0.01*y for y in indexes]
plt.plot(x, cor3[indexes], 'r-')
plt.plot(x, cor1[indexes], 'k--')
plt.plot(x, cor2[indexes], 'b:')
plt.xlabel(r'$r/\sigma$', fontsize=15)
plt.legend(['T = 0.5', 'T = 0.75', 'T = 1.0'], loc='upper right', fontsize=12)

plt.show()
