# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 15:48:34 2023

@author: wattersb
"""

'''
Model the motion of a damped driven non-linear pendulum
'''

import numpy as np
import matplotlib.pyplot as plt

# INITIAL VALUES FOR K PHI AND A
k = 0.5
phi = 0.66667
A = 0.9

# INITIAL VALUE OF dt
dt = 0.01

# INITIAL VALUE OF ITERATION NUMBER
iteration_number = 0

# DEFINE EQUATION 8c FOR NON-LINEAR PENDULUM
def f(theta, omega, t):
    return -np.sin(theta) -k*omega + A*np.cos(phi*t)

###############################
#                             #
#    DIFFERENT VALUES OF A    #
#                             #
###############################

# LIST OF INITIAL CONDITIONS
A_list = np.array([0.90, 1.07, 1.35, 1.47, 1.5])
b_list = np.array([0, 0.1, 0, -0.5, 0])
transient_list = np.array([5000, 5000, 5000, 10000, 100000])
start_graph_list = np.array([10000, 10000, 25000, 30000, 150000])

for x in range(len(A_list)):

    # SET VARIABLES
    theta = 3.0 
    omega = 0.0
    t = 0.0
    A = A_list[x]
    b = b_list[x]
    
    # INITIALISE LISTS
    omega_list = np.array([])
    theta_list = np.array([])
    time_list = np.array([])
    
    # RESET ITERATION NUMBER
    iteration_number = 0
    
    # SET TRANSIENT
    transient = transient_list[x]
    
    for iteration_number in range(start_graph_list[x]):
        # RUNGE-KUTTA METHOD
        k1a = dt * omega
        k1b = dt * f(theta, omega, t)
        k2a = dt * (omega + k1b/2)
        k2b = dt * f(theta + k1a/2, omega + k1b/2, t + dt/2)
        k3a = dt * (omega + k2b/2)
        k3b = dt * f(theta + k2a/2, omega + k2b/2, t + dt/2)
        k4a = dt * (omega + k3b)
        k4b = dt * f(theta + k3a, omega + k3b, t + dt)
        theta = theta + (k1a + 2*k2a + 2*k3a + k4a) / 6
        omega = omega + (k1b + 2*k2b + 2*k3b + k4b) / 6
        
        # INCREASE TIME
        t = t + dt
        
        # CHECK -PI < THETA < PI
        if (theta > np.pi + b) or (theta < -np.pi + b):
            theta -= 2 * np.pi * np.abs(theta) / theta
        
        # APPEND TO ARRAYS IF ITERATION_NUMBER >= TRANSIENT
        if iteration_number >= transient:
            omega_list = np.append(omega_list, omega)
            theta_list = np.append(theta_list, theta)
            time_list = np.append(time_list, t)
        
    # GRAPH
    plt.scatter(omega_list, theta_list, s=0.05)
    plt.ylabel("$\\theta$")
    plt.xlabel("$\omega$")    
    plt.title("Phase portrait of the damped driven non-linear pendulum for A = " + str(A))
    plt.savefig("images/Exercise5"+"."+str(x)+".png", dpi=300)
    plt.show()