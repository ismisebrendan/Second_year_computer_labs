# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 22:11:38 2023

@author: 
"""
import numpy as np
import matplotlib.pyplot as plt

""" 1 """

# f(V) = bV + cV^2
# b = BD
# c = CD^2

B = 1.6*10**(-4) # Ns/m^2
C = 0.25 # Ns^2/m^4

# bV
def b(VD):
    return B * VD

# cV^2
def c(VD):
    return C * VD**2

# f(V)
def f(V):
    return b(V) + c(V)


##########################
#                        #
#   RANGE OF VD VALUES   #
#                        #
##########################


#################
# BOTH RELEVANT #
#################

x = np.arange(0,1*10**(-3),1*10**(-6))

# PLOT bV AND cV AGAINST VD
plt.title("Graph of the contributions to f(V) against VD")
plt.plot(x, b(x), label="bV", color="blue")
plt.plot(x, c(x), label="cV$^2$", color="red")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("bV, cV$^2$")
plt.savefig("images/Exercise1A_together_both.png", dpi=300)
plt.show()

plt.title("Graph of the B against VD")
plt.plot(x, b(x), label="bV", color="blue")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("bV")
plt.savefig("images/Exercise1A_B_both.png", dpi=300)
plt.show()

plt.title("Graph of the C against VD")
plt.plot(x, c(x), label="cV$^2$", color="red")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("cV$^2$")
plt.savefig("images/Exercise1A_C_both.png", dpi=300)
plt.show()

###################
# LINEAR RELEVANT #
###################

x = np.arange(0,1*10**(-5),1*10**(-8))

# PLOT bV AND cV AGAINST VD
plt.title("Graph of the contributions to f(V) against VD")
plt.plot(x, b(x), label="bV", color="blue")
plt.plot(x, c(x), label="cV$^2$", color="red")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("bV, cV$^2$")
plt.savefig("images/Exercise1A_together_linear.png", dpi=300)
plt.show()

plt.title("Graph of the B against VD")
plt.plot(x, b(x), label="bV", color="blue")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("bV")
plt.savefig("images/Exercise1A_B_linear.png", dpi=300)
plt.show()

plt.title("Graph of the C against VD")
plt.plot(x, c(x), label="cV$^2$", color="red")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("cV$^2$")
plt.savefig("images/Exercise1A_C_linear.png", dpi=300)
plt.show()

######################
# QUADRATIC RELEVANT #
######################

x = np.arange(0,1*10**(-1),1*10**(-3))

# PLOT bV AND cV AGAINST VD
plt.title("Graph of the contributions to f(V) against VD")
plt.plot(x, b(x), label="bV", color="blue")
plt.plot(x, c(x), label="cV$^2$", color="red")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("bV, cV$^2$")
plt.savefig("images/Exercise1A_together_quadratic.png", dpi=300)
plt.show()

plt.title("Graph of the B against VD")
plt.plot(x, b(x), label="bV", color="blue")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("bV")
plt.savefig("images/Exercise1A_B_quadratic.png", dpi=300)
plt.show()

plt.title("Graph of the C against VD")
plt.plot(x, c(x), label="cV$^2$", color="red")
plt.legend()
plt.xlabel("V $\\times$ D")
plt.ylabel("cV$^2$")
plt.savefig("images/Exercise1A_C_quadratic.png", dpi=300)
plt.show()


#######################
#                     #
#   DIFFERENT CASES   #
#                     #
#######################


############
# Baseball #
############

D = 0.07 # m
V = 5 # m/s

print("\n Baseball:")
print("Linear term:" + str(b(V)))
print("Quadratic term:" + str(c(V)))
print("Total:" + str(f(V)))

############
# Oil drop #
############

D = 1.5*10**(-6) # m
V = 5*10**(-5) # m/s

print("\n Oil drop:")
print("Linear term:" + str(b(V)))
print("Quadratic term:" + str(c(V)))
print("Total:" + str(f(V)))

############
# Raindrop #
############

D = 10**(-3) # m
V = 1 # m/s

print("\n Rain drop:")
print("Linear term:" + str(b(V)))
print("Quadratic term:" + str(c(V)))
print("Total:" + str(f(V)))

""" 2 """

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

while b >= 10**(-15):
    
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

# RESET V
V = 0 # m/s

# LIST OF TIME VALUES
t = np.arange(0,1,0.0001)

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
m = 10**(-6)

# RESET TIME, SPEED, HEIGHT
t = 0 # s
Y = 5 # m
V = 0 # m/s

# SET LISTS
mass_list = np.array([])
time_list = np.array([])

while m <= 100:
    # RESET TIME, SPEED, HEIGHT
    t = 0 # s
    Y = 5 # m
    V = 0 # m/s
    
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

""" 3 """

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
Vx_list_res, X_list_res, Vy_list_res, Y_list_res, t_list_res = PosResistance(0,0,1,np.pi/4, m)
Vx_list, X_list, Vy_list, Y_list, t_list = Pos(0,0,1,np.pi/4, m)

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
    Vx_list_res, X_list_res, Vy_list_res, Y_list_res, t_list_res = PosResistance(0,0,1,angle[i], m)
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
m = 10e-10 # kg

mass_list = np.array([])

while m <= 10e-4:
    
    # RESET LIST OF X VALUES
    x_list = np.array([])
    
    # CHECK DIFFERENT ANGLES FOR GIVEN MASS
    for i in range(len(angle)):
        Vx_list_res, X_list_res, Vy_list_res, Y_list_res, t_list_res = PosResistance(0,0,1,angle[i], m)
        x_list = np.append(x_list, max(X_list_res))
    
    # MAX ANGLE FOR GIVEN MASS
    maxX = np.where(x_list == max(x_list))
    maxAngle = angle[maxX]
    theta_optimum = np.append(theta_optimum, maxAngle)
    
    # RECORD MASS
    mass_list = np.append(mass_list, m)
    
    # INCREASE THE MASS
    m = m*1.25

plt.plot(np.log10(mass_list), theta_optimum)
plt.scatter(np.log10(mass_list), theta_optimum, s=1, label="Calculated points")
plt.title("Graph of optimum initial angle against mass")
plt.xlabel("log$_{10}$Mass")
plt.ylabel("Optimum nitial angle (rad)")
plt.legend()
plt.savefig("images/Exercise3B_angles_for_masses.png", dpi=300)
plt.show()

""" 4 """

# INITIAL QUANTITIES
B = 1.6*10**(-4) # Ns/m^2
C = 0.25 # Ns^2/m^4
D = 10**(-4) # m
b = B*D
c = C*D**2
rho = 2*10**3 # kg/m^3
m = 4/3 * np.pi * (D/2)**3 * rho # kg
g = 9.81 # m/s^2


########################
#                      #
#   DEFINE FUNCTIONS   #
#                      #
########################


# POSITION WITH AIR RESISTANCE, QUADRATIC DEPENDENCE
def PosResistanceQuadratic(x0, y0, v0, theta, m): # INITIAL X, INITIAL Y, INITIAL SPEED, INITIAL ANGLE, MASS
    
    t = 0
    dt = 0.001
    
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
        dVy = -g*dt - c/m * np.sqrt(Vx**2+Vy**2) * Vx * dt
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
    dt = 0.001
    
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
    dt = 0.001
    
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
theta = np.pi/12 # rad
v = 1 # m/s

i = 1

while m <= 1000:
    while theta <= np.pi/2:
        while v <= 100:
            Vx_list_quad, X_list_quad, Vy_list_quad, Y_list_quad, t_list_quad = PosResistanceQuadratic(0, 0, v, theta, m)
            Vx_list_lin, X_list_lin, Vy_list_lin, Y_list_lin, t_list_lin = PosResistanceLinear(0, 0, v, theta, m)
            Vx_list, X_list, Vy_list, Y_list, t_list = Pos(0, 0, v, theta, m)
            
            plt.plot(X_list_quad, Y_list_quad, label="With air resistance, quadratic dependence")
            plt.plot(X_list_lin, Y_list_lin, label="With air resistance, linear dependence")
            plt.plot(X_list, Y_list, label="Without air resistance")
            plt.legend()
            plt.title("Graph of Y against X for m = " + str(m) + " kg, $\\theta_0$ = " + str(i) + "$\pi/12$ rad, v$_0$ = " + str(v) + " m/s")
            plt.ylabel("Y")
            plt.xlabel("X")
            plt.savefig("images/Exercise4_m"+str(m) + "_theta" + str(i) + "_v" + str(v) + ".png")
            plt.show()
            
            if v < 10:
                v += 1
            else:
                v += 5
         
        i += 1
        theta = theta * i
        v = 1
        
    m = m*100
    theta = np.pi/12
