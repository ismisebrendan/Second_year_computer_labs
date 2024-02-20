# -*- coding: utf-8 -*-
"""
Created on Tue Jan 31 14:55:17 2023

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

# DEFINE THE TOLERANCE OF THE ROOT
tol = 0.0001

# PLOT THE FUNCTION, ITS DERIVATIVE AND X-AXIS
plt.plot(x, f(x))
plt.plot(x, df(x), color='g')
plt.plot(x, 0.0*x, color='black')


# INITITAL GUESS OF ROOT
x1 = 1  # initital estimate of the root
x1 = x1 - f(x1)/df(x1)
plt.scatter(x1, f(x1), color='r')
plt.title("Graph of function showning initial guess of root")
plt.show()


####################
#                  #
#    First Root    #
#                  #
####################

# ITERATE TO FIND ROOT
while np.abs(f(x1)) > tol:
    x1 = x1 - f(x1)/df(x1)

# GRAPH OF THE ROOT
plt.plot(x, f(x))
plt.plot(x, df(x), color='g')
plt.plot(x, 0.0*x, color='black')
plt.scatter(x1, f(x1), color='r')
plt.title("Graph of function showning one root")
plt.savefig("2Graph01.png", dpi=300)
plt.show()

# OUTPUT BEST CALCULATION OF THE ROOT, THE VALUE OF THE FUNCTION HERE, AND ITS ABSOLUTE VALUE
print("x1 =", x1)
print("f(x1) =", f(x1))
print("|f(x1)| =", np.abs(f(x1)))


#####################
#                   #
#    Second Root    #
#                   #
#####################

x2 = -4  # initital estimate of the second root

# ITERATE TO FIND ROOT
while np.abs(f(x2)) > tol:
    x2 = x2 - f(x2)/df(x2)
    
# GRAPH OF THE ROOT
plt.plot(x, f(x))
plt.plot(x, df(x), color='g')
plt.plot(x, 0.0*x, color='black')
plt.scatter(x2, f(x2), color='r')
plt.title("Graph of function showning another root")
plt.savefig("2Graph02.png", dpi=300)
plt.show()

print("x2 =", x2)
print("f(x2) =", f(x2))
print("|f(x2)| =", np.abs(f(x2)))


# GRAPH BOTH ROOTS
plt.plot(x, f(x))
plt.plot(x, df(x), color='g')
plt.plot(x, 0.0*x, color='black')
plt.scatter(x1, f(x1), color='r')
plt.scatter(x2, f(x2), color='r')
plt.title("Graph of function showning both roots")
plt.savefig("2Graph03.png", dpi=300)
plt.show()


########################################################
#                                                      #
#    Dependence of number of steps on the tolerance    #
#                                                      #
########################################################

steps_list = np.array([])   # array of the number of steps 
tol_list = np.array([])     # array of the tolerance values
root_list = np.array([])    # array of the roots calculated

# INITIAL TOLERANCE OF THE ROOT
tol = 0.1    # tolerance of root

while tol >= 10e-15:
    
    # INTITIALISE X VALUE
    x1 = 1
    
    nsteps = 0      # counter for number of steps to find root 
    
    # ITERATE TO FIND THE ROOT
    while np.abs(f(x1)) > tol:
        x1 = x1 - f(x1)/df(x1)
        nsteps += 1
    
    # APPEND TO THE ARRAYS
    steps_list = np.append(steps_list, nsteps)
    tol_list = np.append(tol_list, tol)
    root_list = np.append(root_list, x1)
    
    # DECREASE THE TOLERANCE FOR NEXT TIME
    tol = tol/2

# PLOT OF THE NUMBER OF STEPS TAKEN AGAINST THE LOG BASE 10 OF THE TOLERANCE
plt.scatter(np.log10(tol_list), steps_list)
plt.xlabel("log$_{10}$(tol)")
plt.ylabel("nsteps")
plt.title("Graph of the number of steps required to calculate the root against \n the tolerance in the value of the root")
plt.savefig("2Graph04.png", dpi=300)
plt.show()