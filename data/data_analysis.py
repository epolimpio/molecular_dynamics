import matplotlib as matl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pylab as P

data = np.loadtxt("data/energy.dat", skiprows=5)
correl = np.loadtxt("data/correlation.dat", skiprows=5)
vel_correl = np.loadtxt("data/vel_corr.dat", skiprows=5)

print len(correl)

Ek = data[:,0]
Ep = data[:,1]
Etot = data[:,2]
t = data[:,3]
vcorr = vel_correl[:,0]
modvcorr = vel_correl[:,2]

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

t = range(0,150)
z = range(0,150)
for x in range(0,150):
	t[x] = x*0.004
	z[x] = (1+t[x])**(-1.5)


plt.figure(4)
plt.loglog(t, vcorr[:150], t, modvcorr[:150], t, z)

plt.show()