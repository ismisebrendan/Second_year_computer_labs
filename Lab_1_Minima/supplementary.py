# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 16:29:41 2023

@author: wattersb
"""

import numpy as np
import matplotlib.pyplot as plt


# variables

e = 1.602176634e-19 # C

q1 = 1 * e
q2 = 1 * e
q3 = 1 * e

epsilon0 = 8.85418782e-12 # m-3 kg-1 s4 A2


R = 1e-19       # m


# CHARGES
def V(Q1, Q2, phi):
    return Q1*Q2/(4*np.pi*epsilon0) * 1/(2*R*np.abs(np.sin(phi/2)))

def F(Q1, Q2, phi):
    # F = -dV/dx
    h = 0.00001
    return - ( V(Q1, Q2, phi+h) - V(Q1, Q2, phi) )/h

# phi1 = angle between q1 and q2
phi1 = np.arange(0.0001, 2*np.pi, 0.0001)

plt.plot(phi1, V(q1, q2, phi1))
plt.title("Energy of q1 due to q2 for $\phi_1$")
plt.xlabel("$\phi_1$ (rad)")
plt.ylabel("$V_{12}$ (eV)")
plt.show()

plt.plot(phi1, F(q1, q2, phi1))
plt.title("Force on q1 due to q2 for $\phi_1$")
plt.xlabel("$\phi_1$ (rad)")
plt.ylabel("$F_{12}$ (eV nm$^{-1}$)")
plt.show()
