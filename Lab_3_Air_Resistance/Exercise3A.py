# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 15:37:12 2023

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

A = 2


def Yresistance(t):
    return m*g/b * (-(m/b)*np.exp(-(b*t)/m) - t) + (m**2 * g)/(b**2)


def Xresistance(t):
    return - A * m/b * np.exp(-b/m * t) + A * m/b

def Y(t):
    return -t
    
    

t_list = np.arange(0,1,0.01)
plt.plot(Xresistance(t_list), Yresistance(t_list))
plt.title("Graph of X against Y")
plt.ylabel("Y")
plt.xlabel("X")
plt.show()

plt.plot(t_list, Xresistance(t_list))
plt.title("Graph of X against t")
plt.ylabel("X")
plt.xlabel("T")
plt.show()