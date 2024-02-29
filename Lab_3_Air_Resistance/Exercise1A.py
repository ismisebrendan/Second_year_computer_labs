# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 14:01:25 2023

@author: 
"""

import numpy as np
import matplotlib.pyplot as plt

# f(V) = bV + cV^2
# b = BD
# c = CD^2

B = 1.6*10**(-4) # Ns/m^2
C = 0.25 # Ns^2/m^4

# bV
def b(VD):
    return B * VD

# cV^2
def c(VD):
    return C * VD**2

##########################
#                        #
#   RANGE OF VD VALUES   #
#                        #
##########################

#################
# BOTH RELEVANT #
#################

x = np.arange(0,1*10**(-3),1*10**(-6))

# PLOT bV AND cV AGAINST VD
plt.title("Graph of the contributions to f(V) against VD")
plt.plot(x, b(x), label="bV", color="blue")
plt.plot(x, c(x), label="cV$^2$", color="red")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("bV, cV$^2$")
plt.savefig("images/Exercise1A_together_both.png", dpi=300)
plt.show()

plt.title("Graph of the B against VD")
plt.plot(x, b(x), label="bV", color="blue")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("bV")
plt.savefig("images/Exercise1A_B_both.png", dpi=300)
plt.show()

plt.title("Graph of the C against VD")
plt.plot(x, c(x), label="cV$^2$", color="red")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("cV$^2$")
plt.savefig("images/Exercise1A_C_both.png", dpi=300)
plt.show()

###################
# LINEAR RELEVANT #
###################

x = np.arange(0,1*10**(-5),1*10**(-8))

# PLOT bV AND cV AGAINST VD
plt.title("Graph of the contributions to f(V) against VD")
plt.plot(x, b(x), label="bV", color="blue")
plt.plot(x, c(x), label="cV$^2$", color="red")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("bV, cV$^2$")
plt.savefig("images/Exercise1A_together_linear.png", dpi=300)
plt.show()

plt.title("Graph of the B against VD")
plt.plot(x, b(x), label="bV", color="blue")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("bV")
plt.savefig("images/Exercise1A_B_linear.png", dpi=300)
plt.show()

plt.title("Graph of the C against VD")
plt.plot(x, c(x), label="cV$^2$", color="red")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("cV$^2$")
plt.savefig("images/Exercise1A_C_linear.png", dpi=300)
plt.show()

######################
# QUADRATIC RELEVANT #
######################

x = np.arange(0,1*10**(-1),1*10**(-3))

# PLOT bV AND cV AGAINST VD
plt.title("Graph of the contributions to f(V) against VD")
plt.plot(x, b(x), label="bV", color="blue")
plt.plot(x, c(x), label="cV$^2$", color="red")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("bV, cV$^2$")
plt.savefig("images/Exercise1A_together_quadratic.png", dpi=300)
plt.show()

plt.title("Graph of the B against VD")
plt.plot(x, b(x), label="bV", color="blue")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("bV")
plt.savefig("images/Exercise1A_B_quadratic.png", dpi=300)
plt.show()

plt.title("Graph of the C against VD")
plt.plot(x, c(x), label="cV$^2$", color="red")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("cV$^2$")
plt.savefig("images/Exercise1A_C_quadratic.png", dpi=300)
plt.show()
