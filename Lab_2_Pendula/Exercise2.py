# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 14:54:57 2023

@author: wattersb
"""

'''
Model the motion of a non-linear pendulum
'''

import numpy as np
import matplotlib.pyplot as plt

# INITIAL VALUES FOR K PHI AND A
k = 0.0
phi = 0.66667
A = 0.0

# INITIAL VALUE OF dt
dt = 0.01

# DEFINE EQUATION 8c FOR NON-LINEAR PENDULUM
def f(theta, omega, t):
    return -np.sin(theta) -k*omega + A*np.cos(phi*t)

#############################################
#                                           #
#    DIFFERENT VALUES OF THETA AND OMEGA    #
#                                           #
#############################################

# LIST OF INITIAL CONDITIONS
theta_init = np.array([0.2, 1.0, 3.14, 0.0])
omega_init = np.array([0.0, 0.0, 0.0, 1.0])

for x in range(len(theta_init)):

    # SET VARIABLES
    theta0 = theta_init[x]
    omega0 = omega_init[x]
    theta = theta0
    omega = omega0
    t = 0.0
    nsteps = 0

    # INITIALISE LISTS
    omega_list = np.array([])
    theta_list = np.array([])
    time_list = np.array([])
    
    for i in range(1000):
        # TRAPEZOIDAL RULE
        fndt = f(theta, omega, t) * dt
        theta = theta + omega * dt/2 + (omega+ fndt)*dt/2
        omega = omega + fndt/2 + f(theta, omega + fndt, t + dt) * dt/2
        
        # INCREASE TIME
        t = t + dt
        
        # APPEND TO ARRAYS
        omega_list = np.append(omega_list, omega)
        theta_list = np.append(theta_list, theta)
        time_list = np.append(time_list, t)
    
    # GRAPH V TIME
    plt.plot(time_list, theta_list, label='$\\theta$')
    plt.plot(time_list, omega_list, label='$\omega$')
    plt.legend()
    plt.xlabel("time")
    plt.ylabel("$\\theta$, $\omega$")
    plt.title("Graph of $\\theta$ and $\omega$ v time for $\\theta_0$ =" + str(theta0) + ", $\omega_0$ = " + str(omega0) + " \n for a non-linear pendulum")
    plt.ylim([-np.pi, np.pi])
    plt.savefig("images/Exercise2TogetherTime"+"."+str(x)+".png", dpi=300)
    plt.show()