# -*- coding: utf-8 -*-
"""
Created on Tue Mar 21 16:20:39 2023

@author: wattersb
"""
import numpy as np
import matplotlib.pyplot as plt

# INITIAL QUANTITIES
B = 1.6*10**(-4) # Ns/m^2
D = 10**(-4) # m
b = B*D
rho = 2*10**3 # kg/m^3
m = 4/3 * np.pi * (D/2)**3 * rho # kg
g = 9.81 # m/s^2

def Pos(x0, y0, v0, theta): # INITIAL X, INITIAL Y, INITIAL SPEED, INITIAL ANGLE
    
    t = 0
    dt = 0.0001
    
    X = x0
    Y = y0
    Vx = v0 * np.cos(theta)
    Vy = v0 * np.sin(theta)
    
    # LIST OF POSITIONS, VELOCITIES AND TIMES
    Vx_list = np.array([])
    X_list = np.array([])
    Vy_list = np.array([])
    Y_list = np.array([])
    t_list = np.array([])
    
    while Y >= 0:
        # UPDATE Y STUFF
        dVy = -g*dt
        Vy = Vy + dVy
        Y += Vy*dt
        
        Vy_list = np.append(Vy_list, Vy)
        Y_list = np.append(Y_list, Y)
        
        # UPDATE X STUFF
        X += Vx*dt
        
        Vx_list = np.append(Vx_list, Vx)
        X_list = np.append(X_list, X)
        
        # UPDATE TIME
        t = t + dt
        t_list = np.append(t_list, t)
    
    return Vx_list,X_list,Vy_list,Y_list,t_list

def PosResistance(x0, y0, v0, theta): # INITIAL X, INITIAL Y, INITIAL SPEED, INITIAL ANGLE
    
    t = 0
    dt = 0.0001
    
    X = x0
    Y = y0
    Vx = v0 * np.cos(theta)
    Vy = v0 * np.sin(theta)

    # LIST OF POSITIONS, VELOCITIES AND TIMES
    Vx_list = np.array([])
    X_list = np.array([])
    Vy_list = np.array([])
    Y_list = np.array([])
    t_list = np.array([])
    
    while Y >= 0:
        # UPDATE Y STUFF
        dVy = -g*dt - b/m * Vy * dt
        Vy = Vy + dVy
        Y += Vy*dt
        
        Vy_list = np.append(Vy_list, Vy)
        Y_list = np.append(Y_list, Y)
        
        # UPDATE X STUFF
        dVx = - b/m * Vx * dt
        Vx = Vx + dVx
        X += Vx*dt
        
        Vx_list = np.append(Vx_list, Vx)
        X_list = np.append(X_list, X)
        
        # UPDATE TIME
        t = t + dt
        t_list = np.append(t_list, t)
    
    return Vx_list,X_list,Vy_list,Y_list,t_list

Vx_list_res, X_list_res, Vy_list_res, Y_list_res, t_list_res = PosResistance(0,0,1,np.pi/4)
Vx_list, X_list, Vy_list, Y_list, t_list = Pos(0,0,1,np.pi/4)

plt.plot(X_list_res, Y_list_res, label="Trajectory with air resistance")
plt.plot(X_list, Y_list, label="Trajectory without air resistance")
plt.legend()
plt.title("Graph of Y against X")
plt.ylabel("Y")
plt.xlabel("X")
plt.show()