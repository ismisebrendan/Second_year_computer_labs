# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 15:35:00 2023

@author: wattersb
"""

'''
Model the motion of a damped non-linear pendulum
'''

import numpy as np
import matplotlib.pyplot as plt

# INITIAL VALUES FOR K PHI AND A
k = 0.5
phi = 0.66667
A = 0.0

# INITIAL VALUE OF dt
dt = 0.01

# DEFINE EQUATION 8c FOR NON-LINEAR PENDULUM
def f(theta, omega, t):
    return -np.sin(theta) -k*omega + A*np.cos(phi*t)

# SET VARIABLES
theta = 3.0 
omega = 0.0
t = 0.0

# INITIALISE LISTS
omega_list = np.array([])
theta_list = np.array([])
time_list = np.array([])

for i in range(2000):
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
    
    # APPEND TO ARRAYS
    omega_list = np.append(omega_list, omega)
    theta_list = np.append(theta_list, theta)
    time_list = np.append(time_list, t)

# GRAPHS V TIME
plt.plot(time_list, theta_list)
plt.xlabel("time")
plt.ylabel("$\\theta$")
plt.title("Graph of $\\theta$ v time for $\\theta_0$ = 3.0, $\omega_0$ = 0.0 \n for a damped non-linear pendulum")
plt.ylim([-np.pi, np.pi])
plt.savefig("images/Exercise4Theta.png", dpi=300)
plt.show()

plt.plot(time_list, omega_list)
plt.xlabel("time")
plt.ylabel("$\omega$")
plt.title("Graph of $\omega$ v time for $\\theta_0$ = 3.0, $\omega_0$ = 0.0 \n for a damped non-linear pendulum")
plt.ylim([-np.pi, np.pi])
plt.savefig("images/Exercise4Omega.png", dpi=300)
plt.show()