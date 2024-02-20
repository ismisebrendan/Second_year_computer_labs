# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 16:41:41 2023

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

# GENERATE x POINTS
x = np.arange(-5.0, 2.0, 0.2)

# DEFINE THE TOLERANCE OF THE ROOT
tol = 0.0001

# PLOT THE FUNCTION AND X-AXIS
plt.plot(x, f(x))
plt.plot(x, 0.0*x, color='black')

# GENERATE INITIAL VALUES OF x1 AND x2
x1 = -1      # f(x1) < 0
x3 = 2       # f(x3) > 0

# CHECK INITIAL X VALUES
if f(x1) < 0:
    print("f(x1)<0 -> good")
else:
    print("f(x1)>0 -> not good")
    
if f(x3) > 0:
    print("f(x3)>0 -> good")
else:
    print("f(x3)<0 -> not good")

# DEFINE x2 AS THE MIDPOINT OF x1 AND x3
x2 = 0.5 * (x1 + x3)

plt.scatter(x2, f(x2), color='r') # plot the first guess at the root

# GRAPH OF THE INITIAL GUESS - Part 6
plt.title("Graph of function showning initial midpoint of x1 and x3")
plt.savefig("1.06Graph.png", dpi=300)
plt.show()

####################
#                  #
#    First Root    #
#                  #
####################

# ITERATE TO FIND ROOT
while np.abs(f(x2)) > tol:
    # DEFINE X2 AS THE MIDPOINT
    x2 = 0.5 * (x1 + x3)

    # UPDATE X1 OR X3 TO X2 AND SAVE NEW VALUE
    if f(x2) < 0: # f(x1) < 0 always 
        x1 = x2
    elif f(x2) > 0: # f(x3) > 0 always
        x3 = x2

# OUTPUT BEST CALCULATION OF THE ROOT, THE VALUE OF THE FUNCTION HERE, AND ITS ABSOLUTE VALUE
print("x2 =", x2)
print("f(x2) =", f(x2))
print("|f(x2)| =", np.abs(f(x2)))
    
# GRAPH OF THE ROOT - Part 7  
plt.scatter(x2, f(x2), color='r')   # plot the root
plt.plot(x, f(x))                   # plot the function
plt.plot(x, 0.0*x, color='black')   # plot the x axis
plt.title("Graph of function showning one root")
plt.savefig("1.07Graph.png", dpi=300)
plt.show()

#####################
#                   #
#    Second Root    #
#                   #
#####################

print("---Second Root---")

# GENERATE INITIAL VALUES OF x1 AND x2
x1a = -3      # f(x1) < 0
x3a = -5      # f(x3) > 0

# CHECK INITIAL X VALUES
if f(x1a) < 0:
    print("f(x1)<0 -> good")
else:
    print("f(x1)>0 -> not good")
    
if f(x3a) > 0:
    print("f(x3)>0 -> good")
else:
    print("f(x3)<0 -> not good")

# DEFINE X2 AS THE MIDPOINT
x2a = 0.5 * (x1a + x3a)

# ITERATE TO FIND ROOT
while np.abs(f(x2a)) > tol:
    # DEFINE X2 AS THE MIDPOINT
    x2a = 0.5 * (x1a + x3a)

    # UPDATE X1 OR X3 TO X2 AND OUTPUT NEW VALUES
    if f(x2a) < 0:
        x1a = x2a
    elif f(x2a) > 0:
        x3a = x2a
   
# OUTPUT BEST CALCULATION OF THE ROOT, THE VALUE OF THE FUNCTION HERE, AND ITS ABSOLUTE VALUE
print("x2 =", x2a)
print("f(x2) =", f(x2a))
print("|f(x2)| =", np.abs(f(x2a)))

# GRAPH OF THE ROOT - Part 8  
plt.scatter(x2a, f(x2a), color='r') # plot the root
plt.plot(x, f(x))                   # plot the function
plt.plot(x, 0.0*x, color='black')   # plot the x axis
plt.title("Graph of function showning another root")
plt.show()

# GRAPH OF BOTH ROOTS  
#plot the roots
plt.scatter(x2, f(x2), color='r')   # plot the root
plt.scatter(x2a, f(x2a), color='r') # plot the other root
plt.plot(x, f(x))                   # plot the function
plt.plot(x, 0.0*x, color='black')   # plot the x axis
plt.title("Graph of function showning both roots")
plt.savefig("1.08Graph.png", dpi=300)
plt.show()

########################################################
#                                                      #
#    Dependence of number of steps on the tolerance    #
#                                                      #
########################################################

steps_list = np.array([])   # array of the number of steps 
tol_list = np.array([])     # array of the tolerance values

# INITIAL TOLERANCE OF THE ROOT
tol = 0.1

while tol >= 10e-15:
    
    # INTITIALISE X VALUES
    x1 = -1     
    x3 = 2     
    x2 = 0.5 * (x1 + x3)
    
    nsteps = 0      # counter for number of steps to find root 
    
    # ITERATE TO FIND THE ROOT
    while np.abs(f(x2)) > tol:
        # DEFINE X2 AS THE MIDPOINT
        x2 = 0.5 * (x1 + x3)
    
        # UPDATE X1 OR X3 TO X2 AND OUTPUT NEW VALUES
        if f(x2) < 0:
            x1 = x2
        elif f(x2) > 0:
            x3 = x2    
        nsteps += 1
    
    # APPEND TO THE ARRAYS
    steps_list = np.append(steps_list, nsteps)
    tol_list = np.append(tol_list, tol)
    
    # DECREASE THE TOLERANCE FOR NEXT TIME
    tol = tol/2

# PLOT OF THE NUMBER OF STEPS TAKEN AGAINST THE LOG BASE 10 OF THE TOLERANCE
plt.scatter(np.log10(tol_list), steps_list)
plt.xlabel("log$_{10}$(tol)")
plt.ylabel("nsteps")
plt.title("Graph of the number of steps required to calculate the root \n against the tolerance in the value of the root")
plt.savefig("1.10Graph.png", dpi=300)
plt.show()