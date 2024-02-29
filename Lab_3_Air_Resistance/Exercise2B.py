# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 14:28:51 2023

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
plt.title("Graph of V$_y$ against time")
plt.xlabel("Time (s)")
plt.ylabel("V$_y$ (m/s)")
plt.savefig("images/Exercise2B.png", dpi=300)
plt.show()

############################
#                          #
#   FOR DIFFERENT MASSES   #
#                          #
############################

# INITIAL MASS
m = 10**(-10) # kg

# MAXIMUM TIME
tmax = 10 # s

while m <= 10:
    
    # RESET
    V = 0 # m/s
    t = 0 # s
    V_list = np.array([])
    t_list = np.array([])

    while t < tmax:
        dV = -g*dt - b/m * V * dt
        V = V + dV
        t = t + dt
        
        V_list = np.append(V_list, V)
        t_list = np.append(t_list, t)

    plt.plot(t_list, V_list)
    plt.title("Graph of V$_y$ against time for mass =" + str(m) + "kg")
    plt.xlabel("Time (s)")
    plt.ylabel("V$_y$ (m/s)")
    plt.savefig("images/Exercise2C_mass"+str(m)+".png", dpi=300)
    plt.show()
    
    m = m * 10


#################################
#                               #
#   FOR DIFFERENT RESISTANCES   #
#                               #
#################################


# INITIAL MASS
m = 4/3 * np.pi * (D/2)**3 * rho # kg

# INITIAL RESISTANCE
b = 1.6*10**(-6)

while b >= 10**(-15):
    
    # RESET
    V = 0 # m/s
    t = 0 # s
    V_list = np.array([])
    t_list = np.array([])

    while t < tmax:
        dV = -g*dt - b/m * V * dt
        V = V + dV
        t = t + dt
        
        V_list = np.append(V_list, V)
        t_list = np.append(t_list, t)

    plt.plot(t_list, V_list)
    plt.title("Graph of V$_y$ against time for resistance =" + str(b) + "N s/m")
    plt.xlabel("Time (s)")
    plt.ylabel("V$_y$ (m/s)")
    plt.savefig("images/Exercise2C_resistance"+str(b)+".png", dpi=300)
    plt.show()
    
    b = b / 10




