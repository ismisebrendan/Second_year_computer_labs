# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 21:27:34 2023

@author: wattersb
"""
import numpy as np
import matplotlib.pyplot as plt

# INITIAL QUANTITIES
B = 1.6*10**(-4) # Ns/m^2
C = 0.25 # Ns^2/m^4
D = 10**(-4) # m
b = B*D
c = C*D**2
g = 9.81 # m/s^2


########################
#                      #
#   DEFINE FUNCTIONS   #
#                      #
########################


# POSITION WITH AIR RESISTANCE, QUADRATIC DEPENDENCE
def PosResistanceQuadratic(x0, y0, v0, theta, m): # INITIAL X, INITIAL Y, INITIAL SPEED, INITIAL ANGLE, MASS
    
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
    
    # WHILE ABOVE THE GROUND
    while Y >= 0:
        # UPDATE Y STUFF
        dVy = -g*dt - c/m * np.sqrt(Vx**2+Vy**2) * Vy * dt
        Vy = Vy + dVy
        Y += Vy*dt
        
        Vy_list = np.append(Vy_list, Vy)
        Y_list = np.append(Y_list, Y)
        
        # UPDATE X STUFF
        dVx = -c/m * np.sqrt(Vx**2+Vy**2) * Vx * dt
        Vx = Vx + dVx
        X += Vx*dt
        
        Vx_list = np.append(Vx_list, Vx)
        X_list = np.append(X_list, X)
        
        # UPDATE TIME
        t = t + dt
        t_list = np.append(t_list, t)
    
    return Vx_list, X_list, Vy_list, Y_list, t_list

# POSITION WITH AIR RESISTANCE, LINEAR DEPENDENCE
def PosResistanceLinear(x0, y0, v0, theta, m): # INITIAL X, INITIAL Y, INITIAL SPEED, INITIAL ANGLE, MASS
    
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
    
    # WHILE ABOVE THE GROUND
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
    
    return Vx_list, X_list, Vy_list, Y_list, t_list

# POSITION WITHOUT RESISTANCE
def Pos(x0, y0, v0, theta, m): # INITIAL X, INITIAL Y, INITIAL SPEED, INITIAL ANGLE, MASS
    
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
    
    # WHILE ABOVE THE GROUND
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
    
    return Vx_list, X_list, Vy_list, Y_list, t_list


############################################
#                                          #
#   DIFFERENT MASSES, ANGLES, VELOCITIES   #
#                                          #
############################################

# DEFINE VARIABLES
m = 10**(-10) # kg
theta = 15 # degrees
v = 1 # m/s


while m <= 10**(-2):
    while theta <= 90:
        while v <= 20:
            Vx_list_quad, X_list_quad, Vy_list_quad, Y_list_quad, t_list_quad = PosResistanceQuadratic(0, 0, v, theta * np.pi/180, m)
            Vx_list_lin, X_list_lin, Vy_list_lin, Y_list_lin, t_list_lin = PosResistanceLinear(0, 0, v, theta * np.pi/180, m)
            Vx_list, X_list, Vy_list, Y_list, t_list = Pos(0, 0, v, theta * np.pi/180, m)
            
            plt.plot(X_list_quad, Y_list_quad, label="With air resistance, quadratic dependence")
            plt.plot(X_list_lin, Y_list_lin, label="With air resistance, linear dependence")
            plt.plot(X_list, Y_list, label="Without air resistance")
            plt.legend()
            plt.title("Graph of Y against X for m = " + str(m) + " kg, $\\theta_0$ = " + str(theta) + "$^{\circ}$, v$_0$ = " + str(v) + " m/s")
            plt.ylabel("Y")
            plt.xlabel("X")
            plt.savefig("images/Exercise4_m"+str(m) + "_theta" + str(theta) + "_v" + str(v) + ".png")
            plt.show()
            
            if v < 10:
                v += 1
            else:
                v += 10
    
        theta += 15
        v = 1
        
    m = m*100
    theta = 15