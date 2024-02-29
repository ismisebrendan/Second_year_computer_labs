# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 15:44:12 2023

@author: 
"""

import numpy as np
import matplotlib.pyplot as plt

# f(V) = bV + cV^2
# b = BD
# c = CD^2

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

# MAXIMUM TIME
tmax = 1 # s

########################
#                      #
#   NUMERICAL RESULT   #
#                      #
########################

# LIST OF VELOCITIES AND TIMES
V_list = np.array([])
t_list = np.array([])

while t < tmax:
    dV = -g*dt - b/m * V * dt
    V = V + dV
    t = t + dt
    
    V_list = np.append(V_list, V)
    t_list = np.append(t_list, t)

plt.plot(t_list, V_list)
plt.title("Graph of V$_y$ against time - numerical result")
plt.xlabel("Time (s)")
plt.ylabel("V$_y$ (m/s)")
plt.show()


#########################
#                       #
#   ANALYTICAL RESULT   #
#                       #
#########################

# RESET V
V = 0 # m/s

# LIST OF TIME VALUES
t = np.arange(0,1,0.0001)

def Vy(t):
    return m * g / b * (np.exp(-b*t/m) - 1)

plt.plot(t, Vy(t))
plt.title("Graph of V$_y$ against time - analytical result")
plt.xlabel("Time (s)")
plt.ylabel("V$_y$ (m/s)")
plt.show()


#############
#           #
#   ERROR   #
#           #
#############

plt.plot(t_list, Vy(t_list)-V_list)
plt.title("Graph of the error from the numerical result against time")
plt.xlabel("Time (s)")
plt.ylabel("Error")
plt.show()


