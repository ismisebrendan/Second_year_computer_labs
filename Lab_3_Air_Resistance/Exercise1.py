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

# f(V)
def f(VD):
    return b(VD) + c(VD)


##########################
#                        #
#   RANGE OF VD VALUES   #
#                        #
##########################


#################
# BOTH RELEVANT #
#################

VD = np.arange(0, 10**(-3), 10**(-6))

# PLOT bV AND cV AGAINST VD
plt.title("Graph of the contributions to f(V) against VD")
plt.plot(VD, b(VD), label="bV", color="blue")
plt.plot(VD, c(VD), label="cV$^2$", color="red")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("bV, cV$^2$")
plt.savefig("images/Exercise1A_together_both.png", dpi=300)
plt.show()

plt.title("Graph of the B against VD")
plt.plot(VD, b(VD), label="bV", color="blue")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("bV")
plt.savefig("images/Exercise1A_B_both.png", dpi=300)
plt.show()

plt.title("Graph of the C against VD")
plt.plot(VD, c(VD), label="cV$^2$", color="red")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("cV$^2$")
plt.savefig("images/Exercise1A_C_both.png", dpi=300)
plt.show()

###################
# LINEAR RELEVANT #
###################

VD = np.arange(0, 10**(-5), 10**(-8))

# PLOT bV AND cV AGAINST VD
plt.title("Graph of the contributions to f(V) against VD")
plt.plot(VD, b(VD), label="bV", color="blue")
plt.plot(VD, c(VD), label="cV$^2$", color="red")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("bV, cV$^2$")
plt.savefig("images/Exercise1A_together_linear.png", dpi=300)
plt.show()

plt.title("Graph of the B against VD")
plt.plot(VD, b(VD), label="bV", color="blue")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("bV")
plt.savefig("images/Exercise1A_B_linear.png", dpi=300)
plt.show()

plt.title("Graph of the C against VD")
plt.plot(VD, c(VD), label="cV$^2$", color="red")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("cV$^2$")
plt.savefig("images/Exercise1A_C_linear.png", dpi=300)
plt.show()

######################
# QUADRATIC RELEVANT #
######################

VD = np.arange(0, 10**(-1), 10**(-3))

# PLOT bV AND cV AGAINST VD
plt.title("Graph of the contributions to f(V) against VD")
plt.plot(VD, b(VD), label="bV", color="blue")
plt.plot(VD, c(VD), label="cV$^2$", color="red")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("bV, cV$^2$")
plt.savefig("images/Exercise1A_together_quadratic.png", dpi=300)
plt.show()

plt.title("Graph of the B against VD")
plt.plot(VD, b(VD), label="bV", color="blue")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("bV")
plt.savefig("images/Exercise1A_B_quadratic.png", dpi=300)
plt.show()

plt.title("Graph of the C against VD")
plt.plot(VD, c(VD), label="cV$^2$", color="red")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("cV$^2$")
plt.savefig("images/Exercise1A_C_quadratic.png", dpi=300)
plt.show()


#######################
#                     #
#   DIFFERENT CASES   #
#                     #
#######################


############
# Baseball #
############

D = 0.07 # m
V = 5 # m/s
VD = V*D

print("\n Baseball:")
print("Linear term:" + str(b(VD)))
print("Quadratic term:" + str(c(VD)))
print("Total:" + str(f(VD)))

############
# Oil drop #
############

D = 1.5*10**(-6) # m
V = 5*10**(-5) # m/s
VD = V*D

print("\n Oil drop:")
print("Linear term:" + str(b(VD)))
print("Quadratic term:" + str(c(VD)))
print("Total:" + str(f(VD)))

############
# Raindrop #
############

D = 10**(-3) # m
V = 1 # m/s
VD = V*D

print("\n Rain drop:")
print("Linear term:" + str(b(VD)))
print("Quadratic term:" + str(c(VD)))
print("Total:" + str(f(VD)))
