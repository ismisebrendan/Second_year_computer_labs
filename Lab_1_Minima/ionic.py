# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 15:39:03 2023

@author: wattersb
"""

import numpy as np
import matplotlib.pyplot as plt

# VALUES OF FUNCTION
k = 1.44    # eV nm,   k = e^2/(4*pi*epsilon0)
A = 1090    # eV
p = 0.033   # nm

# POTENTIAL
def V(x):
    return A*np.exp(-x/p) - k/x

# FORCE
def F(x):
    return A/p * np.exp(-x/p) - k * 1/(x**2)

# DERIVATIVE OF FORCE
def dF(x):
    return -A/(p**2) * np.exp(-x/p) + 2*k/(x**3)

###########################
#                         #
#    FIND MINIMUM OF V    #
#                         #
###########################

# FIND ROOT OF F
tol = 10e-10
x1 = 0.2    # initial guess

while np.abs(F(x1)) > tol:
    x1 = x1 - F(x1)/dF(x1)

# PLOT GRAPHS
x = np.arange(0.1, 1.0, 0.001)

plt.plot(x, V(x), label="Potential")
plt.plot(x, F(x), label="Force")
plt.scatter(x1, V(x1), color='black', label="Minimum of potential")
plt.plot(x, np.zeros(x.size), color="black")

plt.title("Plot of length (nm) against energy (eV)")
plt.ylabel("Potential Energy (eV)")
plt.xlabel("Length (nm)")
plt.legend()

plt.ylim([-20,20])
plt.savefig("ionicGraph.png", dpi=300)
plt.show()

# PRINT DATA
print("Minimum of V is at x =", x1, "nm")
print("Value of V at this point is", V(x1), "eV")
print("Value of F at this point is", F(x1), "eV/nm")