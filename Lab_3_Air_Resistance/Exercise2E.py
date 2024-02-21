# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 16:09:18 2023

@author: wattersb
"""

import numpy as np
import matplotlib.pyplot as plt

# f(V) = bV + cV^2
# b = BD

# INITIAL QUANTITIES
B = 1.6*10**(-4) # Ns/m^2
D = 10**(-4) # m
b = B*D
rho = 2*10**3 # kg/m^3
m = 4/3 * np.pi * (D/2)**3 * rho # kg
V = 0 # m/s
t = 0 # s
g = 9.81 # m/s^2
dt = 0.0001

# Y = ANTIDERIVATIVE OF V + c
# c = h0 + gm^2/b^2 if starting from Y = h0 m
# h0 = 5
def Y(t):
    return m*g/b * (-(m/b)*np.exp(-(b*t)/m) - t) + 5 + (m**2 * g)/(b**2)
    

#######################
#                     #
#   TIME TO FALL 5m   #
#                     #
#######################

while Y(t) > 0:
    t = t + dt

print("Time taken to fall from a height of 5 m is " + str(t) + " s")

###########################
#                         #
#   TIME FALLING V MASS   #
#                         #
###########################

# SET MASS
m = 10**(-6)
while Y(t) > 0:
    t = t + dt

print("Time taken to fall from a height of 5 m is " + str(t) + " s")

# SET LISTS
mass_list = np.array([])
time_list = np.array([])


while m <= 100:
    # RESET TIME, SPEED, HEIGHT
    t = 0 # s
    Y = 5
    V = 0
    
    while Y >= 0:
        dV = -g*dt - b/m * V * dt
        V = V + dV
        t = t + dt
        Y += V*dt
                
    
    # APPEND TO ARRAYS
    time_list = np.append(time_list, t)
    mass_list = np.append(mass_list, m)

    # INCREASE MASS
    m = m * 2
    
    

plt.plot(np.log10(mass_list), time_list)
plt.scatter(np.log10(mass_list), time_list)
plt.title("Graph of the time taken to fall 5 m against the mass of a particle")
plt.ylabel("Time (s)")
plt.xlabel("log$_{10}$ mass")
plt.show()



