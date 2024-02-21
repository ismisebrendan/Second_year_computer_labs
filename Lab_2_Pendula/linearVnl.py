# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 15:01:42 2023

@author: wattersb
"""

'''
Compare the motion of a linear and nonlinear pendulum
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
def fNL(theta, omega, t):
    return -np.sin(theta) -k*omega + A*np.cos(phi*t)

# DEFINE EQUATION 8c FOR LINEAR PENDULUM
def fL(theta, omega, t):
    return -theta -k*omega + A*np.cos(phi*t)

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
    t = 0.0
    nsteps = 0

    # SET NONLINEAR AND LINEAR VERSIONS OF VARIABLES
    NL_theta = theta0
    NL_omega = omega0
    
    L_theta = theta0
    L_omega = omega0
    
    # INITIALISE LISTS FOR LINEAR AND NONLINEAR
    NL_omega_list = np.array([])
    NL_theta_list = np.array([])
    
    L_omega_list = np.array([])
    L_theta_list = np.array([])

    time_list = np.array([])
    
    for nsteps in range(1000):
        # TRAPEZOIDAL RULE
        
        # NONLINEAR
        NL_fndt = fNL(NL_theta, NL_omega, t) * dt
        NL_theta = NL_theta + NL_omega * dt/2 + (NL_omega+ NL_fndt)*dt/2
        NL_omega = NL_omega + NL_fndt/2 + fNL(NL_theta, NL_omega + NL_fndt, t + dt) * dt/2
        
        # LINEAR
        L_fndt = fL(L_theta, L_omega, t) * dt
        L_theta = L_theta + L_omega * dt/2 + (L_omega+ L_fndt)*dt/2
        L_omega = L_omega + L_fndt/2 + fL(L_theta, L_omega + L_fndt, t + dt) * dt/2
        
        # INCREASE TIME
        t = t + dt
        
        # APPEND TO ARRAYS
        NL_omega_list = np.append(NL_omega_list, NL_omega)
        NL_theta_list = np.append(NL_theta_list, NL_theta)
        
        L_omega_list = np.append(L_omega_list, L_omega)
        L_theta_list = np.append(L_theta_list, L_theta)
        
        time_list = np.append(time_list, t)
    
    # GRAPHs V TIME
    plt.plot(time_list, NL_theta_list, label='$\\theta$ Non linear')
    plt.plot(time_list, L_theta_list, label='$\\theta$ Linear')
    plt.legend()
    plt.xlabel("time")
    plt.ylabel("$\\theta$")
    plt.title("Graph of $\\theta$ v time for $\\theta_0$ =" + str(theta0) + ", $\omega_0$ = " + str(omega0) + " \n for non-linear and linear pendula")
    plt.ylim([-np.pi, np.pi])
    plt.savefig("images/ComparisonTheta"+"."+str(i)+".png", dpi=300)
    plt.show()

    plt.plot(time_list, NL_omega_list, label='$\omega$ Non linear')
    plt.plot(time_list, L_omega_list, label='$\omega$ Linear')
    plt.legend()
    plt.xlabel("time")
    plt.ylabel("$\\theta$, $\omega$")
    plt.title("Graph of $\omega$ v time for $\\theta_0$ =" + str(theta0) + ", $\omega_0$ = " + str(omega0) + " \n for non-linear and linear pendula")
    plt.ylim([-np.pi, np.pi])
    plt.savefig("images/ComparisonOmega"+"."+str(i)+".png", dpi=300)
    plt.show()