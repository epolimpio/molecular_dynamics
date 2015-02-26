import matplotlib as matl
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import pylab as P

data = np.loadtxt("energy006.dat", skiprows=5)
correl = np.loadtxt("correlation006.dat", skiprows=5)
vel_correl = np.loadtxt("vel_corrGas.dat", skiprows=5)
p = np.loadtxt("pressure006.dat", skiprows=5)

print 1+1/(3*864*0.5)*np.mean(p)+16*np.pi*1.2/(3*0.5)*(2/(3*3.3**9) - 1/(3.3**3))

Ek = data[:,0]
Ep = data[:,1]
Etot = data[:,2]
t = data[:,3]
vcorr = vel_correl[:,0]
modvcorr = vel_correl[:,2]

print 1/(1/(1.5*864)-np.var(Ek)/np.mean(Ek)**2)

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


plt.figure(4)
t = range(200,350)
for x in range(0,150):
	t[x] = t[x]*0.004
z = np.polyfit(np.log(t),np.log(modvcorr[200:350]),1)

print z
plt.loglog(range(100,350), modvcorr[100:350])

plt.show()