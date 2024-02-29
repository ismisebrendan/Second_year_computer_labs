# -*- coding: utf-8 -*-
"""
Created on Tue Mar 14 14:28:51 2023

@author: 
"""
import numpy as np
import matplotlib.pyplot as plt

# f(V) = bV + cV^2
# b = BD
# c = CD^2

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

# MAXIMUM TIME
tmax = 1 # s

# LIST OF VELOCITIES AND TIMES
V_list = np.array([])
t_list = np.array([])

while t < tmax:
    dV = -g*dt - b/m * V * dt
    V = V + dV
    t = t + dt
    
    V_list = np.append(V_list, V)
    t_list = np.append(t_list, t)

plt.plot(t_list, V_list)
plt.title("Graph of V$_y$ against time")
plt.xlabel("Time (s)")
plt.ylabel("V$_y$ (m/s)")
plt.savefig("images/Exercise2B.png", dpi=300)
plt.show()


############################
#                          #
#   FOR DIFFERENT MASSES   #
#                          #
############################


# INITIAL MASS
m = 10**(-10) # kg

# MAXIMUM TIME
tmax = 10 # s

while m <= 10:
    
    # RESET
    V = 0 # m/s
    t = 0 # s
    V_list = np.array([])
    t_list = np.array([])

    while t < tmax:
        dV = -g*dt - b/m * V * dt
        V = V + dV
        t = t + dt
        
        V_list = np.append(V_list, V)
        t_list = np.append(t_list, t)

    plt.plot(t_list, V_list)
    plt.title("Graph of V$_y$ against time for mass =" + str(m) + "kg")
    plt.xlabel("Time (s)")
    plt.ylabel("V$_y$ (m/s)")
    plt.savefig("images/Exercise2C_mass"+str(m)+".png", dpi=300)
    plt.show()
    
    m = m * 10


#################################
#                               #
#   FOR DIFFERENT RESISTANCES   #
#                               #
#################################


# INITIAL MASS
m = 4/3 * np.pi * (D/2)**3 * rho # kg

# INITIAL RESISTANCE
b = 1.6*10**(-6)

while b >= 1.6*10**(-15):
    
    # RESET
    V = 0 # m/s
    t = 0 # s
    V_list = np.array([])
    t_list = np.array([])

    while t < tmax:
        dV = -g*dt - b/m * V * dt
        V = V + dV
        t = t + dt
        
        V_list = np.append(V_list, V)
        t_list = np.append(t_list, t)

    plt.plot(t_list, V_list)
    plt.title("Graph of V$_y$ against time for resistance =" + str(b) + "N s/m")
    plt.xlabel("Time (s)")
    plt.ylabel("V$_y$ (m/s)")
    plt.savefig("images/Exercise2C_resistance"+str(b)+".png", dpi=300)
    plt.show()
    
    b = b / 10


################################################
#                                              #
#   COMPARE NUMERICAL AND ANALYTICAL RESULTS   #
#                                              #
################################################


####################
# NUMERICAL RESULT #
####################

# LIST OF VELOCITIES AND TIMES
V_list = np.array([])
t_list = np.array([])

tmax = 1 # s

# RESET V, TIME
V = 0 # m/s
t = 0 # s
b = B*D
m = 4/3 * np.pi * (D/2)**3 * rho # kg

while t < tmax:
    dV = -g*dt - b/m * V * dt
    V = V + dV
    t = t + dt
    
    V_list = np.append(V_list, V)
    t_list = np.append(t_list, t)

plt.plot(t_list, V_list)
plt.title("Graph of V$_y$ against time - numerical result")
plt.xlabel("Time (s)")
plt.ylabel("V$_y$ (m/s)")
plt.savefig("images/Exercise2D_numerical.png", dpi=300)
plt.show()

#####################
# ANALYTICAL RESULT #
#####################

# RESET VARIABLES
V = 0 # m/s
b = B*D
m = 4/3 * np.pi * (D/2)**3 * rho # kg

# LIST OF TIME VALUES
t = np.arange(0, 1, 0.0001)

def Vy(t):
    return m * g / b * (np.exp(-b*t/m) - 1)

plt.plot(t, Vy(t))
plt.title("Graph of V$_y$ against time - analytical result")
plt.xlabel("Time (s)")
plt.ylabel("V$_y$ (m/s)")
plt.savefig("images/Exercise2D_analytical.png", dpi=300)
plt.show()

#########
# ERROR #
#########

plt.plot(t_list, Vy(t_list)-V_list)
plt.title("Graph of the error from the numerical result against time")
plt.xlabel("Time (s)")
plt.ylabel("Error")
plt.savefig("images/Exercise2D_error.png", dpi=300)
plt.show()


###########################
#                         #
#   TIME FALLING V MASS   #
#                         #
###########################

# SET MASS
m = 10**(-8)

# SET LISTS
mass_list = np.array([])
time_list = np.array([])

while m <= 100:
    # RESET TIME, SPEED, HEIGHT
    t = 0 # s
    Y = 5
    V = 0
    
    while Y >= 0:
        dV = -g*dt - b/m * V * dt
        V = V + dV
        t = t + dt
        Y += V*dt
                
    # APPEND TO ARRAYS
    time_list = np.append(time_list, t)
    mass_list = np.append(mass_list, m)

    # INCREASE MASS
    m = m * 2
    
plt.plot(np.log10(mass_list), time_list)
plt.scatter(np.log10(mass_list), time_list)
plt.title("Graph of the time taken to fall 5 m against the mass of a particle")
plt.ylabel("Time (s)")
plt.xlabel("log$_{10}$ mass")
plt.savefig("images/Exercise2E.png", dpi=300)
plt.show()
