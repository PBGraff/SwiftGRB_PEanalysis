import numpy as np
from scipy.integrate import quad, romberg

MPC_IN_M = 3.08567758e22
GPC_IN_M = (1e3 * MPC_IN_M)
C_LIGHT = 2.99792458e8
HUBBLE0 = (1e5 * 0.71 / MPC_IN_M)
DH = (C_LIGHT / HUBBLE0 / GPC_IN_M)
DH3 = (DH * DH * DH)

Omega_m = 0.274
Omega_lambda = 0.726

z1 = 3.6
n0 = 0.84
n1 = 2.07
n2 = -0.7

def E(z):
	return np.sqrt(Omega_m * np.power((1.0 + z), 3.0) + Omega_lambda)

def Einv(z):
	return 1.0 / E(z)

def Einvint(z):
	return romberg(Einv, 0, z)

def R(z):
	if z<z1:
		return n0 * np.power(1.0+z, n1)
	else:
		return n0 * np.power(1.0+z1, n1-n2) * np.power(1.0+z, n2)

def Rdz(z):
	return R(z) / (1.0 + z) * DH3 / E(z) * np.power(Einvint(z), 2.0)

Number = 4.0 * np.pi * romberg(Rdz, 0.0, 10.0)

print Number
