import matplotlib as matl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pylab as P

data = np.loadtxt("data/energy02.dat", skiprows=3)
correl = np.loadtxt("data/correlation02.dat", skiprows=3)

print len(correl)

Ek = data[:,0]
Ep = data[:,1]
Etot = data[:,2]
t = data[:,3]

n = range(0,len(Ek))

plt.figure(1)
plt.subplot(311)
plt.plot(n, Ek, 'b', label = 'Kinetic')

plt.subplot(312)
plt.plot(n, Ep, 'r', label = 'Potential')

plt.subplot(313)
plt.plot(n, Etot, 'g', label = 'Total')

plt.figure(2)
plt.plot(n, t, 'g', label = 'Temperature')

plt.figure(3)
plt.plot(range(0,len(correl)), correl, 'g', label = 'Temperature')

plt.show()