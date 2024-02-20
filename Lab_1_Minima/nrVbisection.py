# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 15:22:21 2023

@author: wattersb
"""

import numpy as np
import matplotlib.pyplot as plt

# COEFFICIENTS OF THE PARABOLA
a = 2
b = 5
c = -10

# DEFINE PARABOLIC FUNCTION
def f(x):
    return a*x**2 + b*x + c

# DEFINE DERIVATIVE
def df(x):       # x position, h
    return 2*a*x + b

# GENERATE x POINTS
x = np.arange(-5.0, 2.0, 0.2)


# INITIALISE ARRAYS
steps_list_nr = np.array([])          # array of the number of steps for Newton-Raphson Method
steps_list_bisection = np.array([])   # array of the number of steps for Bisection Method
tol_list = np.array([])               # array of the tolerance values

# INITIAL TOLERANCE OF THE ROOT
tol = 0.1

while tol >= 10e-15:
    
    ##########################
    #                        #
    #    Bisection Method    #
    #                        #
    ##########################
    
    # INTITIALISE X VALUES
    x1 = -1     
    x3 = 2 
    x2 = 0.5 * (x1 + x3)
    
    nsteps_bisecton = 0      # counter for number of steps to find root 
    
    # ITERATE TO FIND THE ROOT
    while np.abs(f(x2)) > tol:
        # DEFINE X2 AS THE MIDPOINT
        x2 = 0.5 * (x1 + x3)
    
        # UPDATE X1 OR X3 TO X2 AND OUTPUT NEW VALUES
        if f(x2) < 0:
            x1 = x2
        elif f(x2) > 0:
            x3 = x2    
        nsteps_bisecton += 1
    
    ###############################
    #                             #
    #    Newton-Raphson Method    #
    #                             #
    ###############################
    # INTITIALISE X VALUE
    x1 = 1
    
    nsteps_nr = 0      # counter for number of steps to find root 
    
    # ITERATE TO FIND THE ROOT
    while np.abs(f(x1)) > tol:
        x1 = x1 - f(x1)/df(x1)
        nsteps_nr += 1
    
    # APPEND TO THE ARRAYS
    steps_list_bisection = np.append(steps_list_bisection, nsteps_bisecton)
    steps_list_nr = np.append(steps_list_nr, nsteps_nr)
    tol_list = np.append(tol_list, tol)
    
    # DECREASE THE TOLERANCE FOR NEXT TIME
    tol = tol/2

# GRAPH THE TWO METHODS AGAINST EACH OTHER
plt.scatter(np.log10(tol_list), steps_list_bisection, color='b', label="Bisection method")
plt.scatter(np.log10(tol_list), steps_list_nr, color='g', label="Newton-Raphson method")
plt.legend()

plt.xlabel("log$_{10}$(tol)")
plt.ylabel("nsteps")
plt.title("Graph of the number of steps required to calculate the root against \n the tolerance in the value of the root for each root")
plt.savefig("nrVbisection.png", dpi=300)
plt.show()