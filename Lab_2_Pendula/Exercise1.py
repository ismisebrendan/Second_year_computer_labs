# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 13:59:32 2023

@author: wattersb
"""

'''
Model the motion of a linear pendulum
'''

import numpy as np
import matplotlib.pyplot as plt

# INITIAL VALUES FOR K PHI AND A
k = 0.0
phi = 0.66667
A = 0.0

# INITIAL VALUES OF VARIABLES
theta = 0.2
omega = 0.0
t = 0.0
dt = 0.01
nsteps = 0

# DEFINE EQUATION 8c FOR LINEAR PENDULUM
def f(theta, omega, t):
    return -theta -k*omega + A*np.cos(phi*t)

# CHECK FUNCTION BY PRINTING f FOR DIFFERENT THETA, OMEGA, T VALUES -> EFFECTIVLY PRINT -THETA
print(f(theta, omega, t))
print(f(0.1, 0.56, 2))
print(f(theta, 20, t))

#############################################
#                                           #
#    DIFFERENT VALUES OF THETA AND OMEGA    #
#                                           #
#############################################

# LIST OF INITIAL CONDITIONS
theta_init = np.array([0.2, 1.0, 3.14, 0.0])
omega_init = np.array([0.0, 0.0, 0.0, 1.0])

for i in range(len(theta_init)):

    # SET VARIABLES
    theta0 = theta_init[i]
    omega0 = omega_init[i]
    theta = theta0
    omega = omega0
    t = 0.0
    nsteps = 0
    
    # INITIALISE LISTS
    omega_list = np.array([])
    theta_list = np.array([])
    step_list = np.array([])
    time_list = np.array([])
    
    # TRAPEZOIDAL RULE
    for nsteps in range(1000):
        fndt = f(theta, omega, t) * dt
        theta = theta + omega * dt/2 + (omega+ fndt)*dt/2
        omega = omega + fndt/2 + f(theta, omega + fndt, t + dt) * dt/2
        
        # INCREASE TIME
        t = t + dt
        
        # APPEND TO ARRAYS
        omega_list = np.append(omega_list, omega)
        theta_list = np.append(theta_list, theta)
        step_list = np.append(step_list, nsteps)
        time_list = np.append(time_list, t)
    
    # SEPARATE GRAPHS
    plt.plot(time_list, theta_list)
    plt.xlabel("time")
    plt.ylabel("$\\theta$")
    plt.title("Graph of $\\theta$ v time for $\\theta_0$ =" + str(theta0) + ", $\omega_0$ = " + str(omega0) + " \n for a linear pendulum")
    plt.ylim([-np.pi, np.pi])
    plt.savefig("images/ThetaVtime"+str(theta0)+"."+str(i)+".png", dpi=300)
    plt.show()
    
    plt.plot(time_list, omega_list)
    plt.xlabel("time")
    plt.ylabel("$\omega$")
    plt.title("Graph of $\omega$ v time for $\\theta_0$ =" + str(theta0) + ", $\omega_0$ = " + str(omega0) + " \n for a linear pendulum")
    plt.ylim([-np.pi, np.pi])
    plt.savefig("images/OmegaVtime"+str(omega0)+"."+str(i)+".png", dpi=300)
    plt.show()
    
    # GRAPHS TOGETHER
    plt.plot(step_list, theta_list, label='$\\theta$')
    plt.plot(step_list, omega_list, label='$\omega$')
    plt.legend()
    plt.xlabel("nsteps")
    plt.ylabel("$\\theta$, $\omega$")
    plt.title("Graph of $\\theta$ and $\omega$ v nsteps for $\\theta_0$ =" + str(theta0) + ", $\omega_0$ = " + str(omega0) + " \n for a linear pendulum")
    plt.axis([0, 500, -np.pi, np.pi])
    plt.savefig("images/Exercise1TogetherSteps"+"."+str(i)+".png", dpi=300)
    plt.show()
    
    # GRAPH V TIME
    plt.plot(time_list, theta_list, label='$\\theta$')
    plt.plot(time_list, omega_list, label='$\omega$')
    plt.legend()
    plt.xlabel("time")
    plt.ylabel("$\\theta$, $\omega$")
    plt.title("Graph of $\\theta$ and $\omega$ v time for $\\theta_0$ =" + str(theta0) + ", $\omega_0$ = " + str(omega0) + " \n for a linear pendulum")
    plt.ylim([-np.pi, np.pi])
    plt.savefig("images/Exercise1TogetherTime"+"."+str(i)+".png", dpi=300)
    plt.show()