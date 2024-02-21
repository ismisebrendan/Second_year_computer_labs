# -*- coding: utf-8 -*-
"""
Created on Tue Feb 21 15:19:59 2023

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

# SET VARIABLES
theta = 3.14
omega = 0.0
t = 0.0

RK_theta = theta
RK_omega = omega

TR_theta = theta
TR_omega = omega

# INITIALISE LISTS
RK_theta_list = np.array([])
TR_theta_list = np.array([])
time_list = np.array([])

for i in range(1000):
    # RUNGE-KUTTA
    k1a = dt * RK_omega
    k1b = dt * f(RK_theta, RK_omega, t)
    k2a = dt * (RK_omega + k1b/2)
    k2b = dt * f(RK_theta + k1a/2, RK_omega + k1b/2, t + dt/2)
    k3a = dt * (RK_omega + k2b/2)
    k3b = dt * f(RK_theta + k2a/2, RK_omega + k2b/2, t + dt/2)
    k4a = dt * (RK_omega + k3b)
    k4b = dt * f(RK_theta + k3a, RK_omega + k3b, t + dt)
    RK_theta = RK_theta + (k1a + 2*k2a + 2*k3a + k4a) / 6
    RK_omega = RK_omega + (k1b + 2*k2b + 2*k3b + k4b) / 6
    
    # TRAPEZOIDAL RULE
    fndt = f(TR_theta, TR_omega, t) * dt
    TR_theta = TR_theta + TR_omega * dt/2 + (TR_omega+ fndt)*dt/2
    TR_omega = TR_omega + fndt/2 + f(TR_theta, TR_omega + fndt, t + dt) * dt/2
    
    # INCREASE TIME
    t = t + dt
    
    # APPEND TO ARRAYS
    RK_theta_list = np.append(RK_theta_list, RK_theta)
    TR_theta_list = np.append(TR_theta_list, TR_theta)
    time_list = np.append(time_list, t)

# GRAPH V TIME
plt.plot(time_list, RK_theta_list, label='$\\theta$ for Runge-Kutta')
plt.plot(time_list, TR_theta_list, label='$\\theta$ for Trapezoidal Rule')
plt.legend()
plt.xlabel("time")
plt.ylabel("$\\theta$")
plt.title("Graph of $\\theta$ and $\omega$ v time for $\\theta_0$ = 3.14, $\omega_0$ = 0.0 for a non-linear \n pendulum calculated using the RK Method and Trapezoidal Rule")
plt.ylim([-np.pi, np.pi])
plt.savefig("images/Exercise3.png", dpi=300)
plt.show()