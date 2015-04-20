import numpy as np
from scipy.integrate import quad
import matplotlib.pyplot as plt

def Redshift(n0, n1, n2, z1=3.6, z=np.linspace(0,10,num=1001)):
	Rlow = np.power((1.0 + z), n1)
	Rhigh = np.power((1.0 + z), n2)
	rbrk = np.power((1.0 + z1), n1 - n2)
	R = Rlow * (z <= z1) + rbrk * Rhigh * (z > z1)
	R *= n0 / R[0]
	return z, R

z, R = Redshift(0.84, 2.07, -0.7)

plt.plot(z,R,'-k')
plt.xlabel(r'$z$')
plt.ylabel(r'$\mathcal{R}(z)$')
plt.grid()
#plt.gca().set_yscale('log')
plt.show()



#### This computes E(z) and int_0^z dz'/E(z') and saves to file

def Efunc(z):
	Omega_m = 0.274
	Omega_lambda = 0.726
	E = np.sqrt(Omega_m * np.power((1 + z), 3) + Omega_lambda)
	return E

def Efuncinv(z):
	return 1.0 / Efunc(z)

z = np.linspace(0,10,num=1001)
dz = z[1] - z[0]

E = Efunc(z)

Eics = np.zeros(E.shape)
for i in range(len(Eics)):
	Eics[i] = (quad(Efuncinv, 0, z[i])[0])**2.0

#Eics = np.square(np.cumsum(1.0 / E) * dz)
#Eics[1:] = Eics[:-1]
#Eics[0] = 0

Eall = Eics / E;

z = z.reshape(z.shape[0],1)
E = E.reshape(E.shape[0],1)
Eics = Eics.reshape(Eics.shape[0],1)
Eall = Eall.reshape(Eall.shape[0],1)

d = np.concatenate((z,E,Eics,Eall),axis=1)

np.savetxt('support_data/splines_Ez.txt',d,fmt='%0.9lf')

z2, R = Redshift(0.84, 2.07, -0.7, z=z)

Rp = R / (1+z) * Eall

plt.plot(z,Rp,'-k')
plt.plot(z,R/(1+z),'--b')
plt.plot(z,Eall,'-.r')
#plt.plot(z,np.cumsum(Eall),'-g')
plt.xlabel(r'$z$')
plt.grid()
plt.show()
