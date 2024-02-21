# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 14:15:54 2023

@author: wattersb
"""

# f(V) = bV + cV^2
# b = BD
# c = CD^2

B = 1.6*10**(-4) # Ns/m^2
C = 0.25 # Ns^2/m^4

# bV
def b(V):
    return B * V * D

# cV^2
def c(V):
    return C * V**2 * D**2

# f(V)
def f(V):
    return b(V) + c(V)

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

print("\n Baseball:")
print("Linear term:" + str(b(V)))
print("Quadratic term:" + str(c(V)))
print("Total:" + str(f(V)))

############
# Oil drop #
############

D = 1.5*10**(-6) # m
V = 5*10**(-5) # m/s

print("\n Oil drop:")
print("Linear term:" + str(b(V)))
print("Quadratic term:" + str(c(V)))
print("Total:" + str(f(V)))

############
# Raindrop #
############

D = 10**(-3) # m
V = 1 # m/s

print("\n Rain drop:")
print("Linear term:" + str(b(V)))
print("Quadratic term:" + str(c(V)))
print("Total:" + str(f(V)))

"""
Ideal forms of f(V) for these situations:
                                                        DV
    Baseball - Quadratic term is sufficient, c(V)       0.35 m^2/s
    Oil drop - Linear term is sufficient, b(V)          7.5e-11 m^2/s
    Raindrop - Total is needed, f(V)                    10^-3 m^2/s
"""