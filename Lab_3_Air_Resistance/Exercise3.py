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


###############################
#                             #
#   TRAJECTORIES OF OBJECTS   #
#                             #
###############################


####################
# DEFINE FUNCTIONS #
####################

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

# POSITION WITH AIR RESISTANCE
def PosResistance(x0, y0, v0, theta, m): # INITIAL X, INITIAL Y, INITIAL SPEED, INITIAL ANGLE, MASS
    
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

# GET THE VALUES
Vx_list_res, X_list_res, Vy_list_res, Y_list_res, t_list_res = PosResistance(0, 0, 1, np.pi/4, m)
Vx_list, X_list, Vy_list, Y_list, t_list = Pos(0, 0, 1, np.pi/4, m)

plt.plot(X_list_res, Y_list_res, label="Trajectory with air resistance")
plt.plot(X_list, Y_list, label="Trajectory without air resistance")
plt.legend()
plt.title("Graph of Y against X")
plt.ylabel("Y")
plt.xlabel("X")
plt.savefig("images/Exercise3A.png", dpi=300)
plt.show()


###############################################
#                                             #
#   TEST DIFFERENT ANGLES FOR CONSTANT MASS   #
#                                             #
###############################################


# SET THE MASS
m = 4/3 * np.pi * (D/2)**3 * rho # kg

# SET LISTS OF ANGLES AND HORIZONTAL DISTANCES
angle = np.arange(0,np.pi/2,0.01)
x_list = np.array([])

for i in range(len(angle)):
    Vx_list_res, X_list_res, Vy_list_res, Y_list_res, t_list_res = PosResistance(0, 0, 1, angle[i], m)
    x_list = np.append(x_list, max(X_list_res))

plt.plot(angle, x_list)
plt.title("Graph of horizontal displacement against initial angle")
plt.ylabel("Horizontal displacement (m)")
plt.xlabel("Initial angle (rad)")

# FIND THE ANGLE OF MAXIMUM DISPLACEMENT AND PRINT AND PLOT IT
maxX = np.where(x_list == max(x_list))
maxAngle = angle[maxX]
print("The maximum displacement occurs with an initial angle of " + str(maxAngle) + " radians")
plt.scatter(maxAngle, x_list[maxX], color='black', label="Angle of Maximum Displacement")
plt.legend()
plt.savefig("images/Exercise3B_constant_mass.png", dpi=300)
plt.show()


###############################################
#                                             #
#   TEST DIFFERENT ANGLES FOR CHANGING MASS   #
#                                             #
###############################################


# LIST OF OPTIMUM ANGLES
theta_optimum = np.array([])

# INITIAL MASS
m = 10**(-9) # kg

mass_list = np.array([])

while m <= 10**(-3):
    
    # RESET LIST OF X VALUES
    x_list = np.array([])
    
    # CHECK DIFFERENT ANGLES FOR GIVEN MASS
    for i in range(len(angle)):
        Vx_list_res, X_list_res, Vy_list_res, Y_list_res, t_list_res = PosResistance(0, 0, 1,angle[i], m)
        x_list = np.append(x_list, max(X_list_res))
    
    # MAX ANGLE FOR GIVEN MASS
    maxX = np.where(x_list == max(x_list))
    maxAngle = angle[maxX]
    theta_optimum = np.append(theta_optimum, maxAngle)
    
    # RECORD MASS
    mass_list = np.append(mass_list, m)
    
    # INCREASE THE MASS
    m = m*2

plt.plot(np.log10(mass_list), theta_optimum)
plt.plot(np.log10(mass_list), np.pi/4 * np.ones(len(mass_list)), color="black", label="$\\theta_0$ = $\\frac{\pi}{4}$")
plt.scatter(np.log10(mass_list), theta_optimum, s=5, label="Calculated points")
plt.title("Graph of optimum initial angle against mass")
plt.xlabel("log$_{10}$Mass")
plt.ylabel("Optimum initial angle (rad)")
plt.legend()
plt.savefig("images/Exercise3B_angles_for_masses.png", dpi=300)
plt.show()