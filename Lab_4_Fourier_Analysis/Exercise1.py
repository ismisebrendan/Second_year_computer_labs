# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 14:05:30 2023

@author: 
"""
import numpy as np
import matplotlib.pyplot as plt

#####################################
#                                   #
#   SIMPSONS RULE FOR INTEGRATION   #
#                                   #
#####################################


def Simpson(f, a, b, n): # FUNCTION, LOWER LIMIT, UPPER LIMIT, NUMBER OF STEPS
    n = int(n)*2 # ENSURE ONLY EVEN NUMBERS OF STEPS ARE USED
    h = (b-a) / n
    
    # THE SUMMATION TERMS
    sigma1 = 0
    sigma2 = 0
    
    for j in range(1, int(n/2)): # LOOP OVER [1,n/2 - 1]
        x2j = a + 2*j*h
        sigma1 += f(x2j)

    for j in range(1, int(n/2 + 1)): # LOOP OVER [1,n/2]
        xj = a + (2*j-1)*h
        sigma2 += f(xj)
    
    return h/3 * (f(a) + 2*sigma1 + 4*sigma2 + f(b))

print("--- Testing the Simpson's Rule code ---")

print("Integral evaluated using Simpson's Rule "+str(Simpson(np.exp, 0, 1, 4)))

print("Difference between Simpson's Rule and actual value "+str(Simpson(np.exp, 0, 1, 4) - (np.exp(1) - np.exp(0))))


#############################################
#                                           #
#   DEFINE FUNCTIONS TO FIND COEFFICIENTS   #
#                                           #
#############################################


def a0(f, omega, n): # NATURAL FREQUENCY, FUNCTION, NUMBER OF STEPS/2
    T = 2*np.pi/omega
    return 1/T * Simpson(f, 0, T, n)

def ak(f, omega, n, k): # NATURAL FREQUENCY, FUNCTION, NUMBER OF STEPS/2, k STEPS
    T = 2*np.pi/omega
    
    # EMPTY ARRAY OF ak VALUES
    ak = np.array([])
    
    # THE FUNCTION TIMES THE COSINE TERM
    def fcos(t):
        return f(t)*np.cos(t*omega*i)
    
    # CALCULATE THE ak TERMS
    for i in range(1, k+1):
        a = Simpson(fcos, 0, T, n)
        
        ak = np.append(ak, a)
            
    ak = 2/T * ak
    return ak

def bk(f, omega, n, k): # NATURAL FREQUENCY, FUNCTION, NUMBER OF STEPS/2, k STEPS
    T = 2*np.pi/omega
    
    # EMPTY ARRAY OF bk VALUES
    bk = np.array([])
    
    # THE FUNCTION TIMES THE SINE TERM
    def fsin(t):
        return f(t)*np.sin(t*omega*i)
    
    # CALCULATE THE bk TERMS
    for i in range(1, k+1):
        b = Simpson(fsin, 0, T, n)
        
        bk = np.append(bk, b)

    bk = 2/T * bk
    return bk

def fourierSeries(f, omega, n, k, t):
    
    # CALCULATE THE COEFFICIENTS
    azero = a0(f, omega, n)
    a = ak(f, omega, n, k)
    b = bk(f, omega, n, k)
    
    # CALCULATE THE SERIES
    fourier = 0
    for m in range(1, len(a)+1):
        fourier += a[m-1]*np.cos(m*t) + b[m-1]*np.sin(m*t)
    fourier += azero
    
    return azero, a, b, fourier


###############################
#                             #
#   FIND THE FOURIER SERIES   #
#                             #
###############################


# FUNCTION TO GRAPH THE FUNCTION, SERIES AND COEFFICIENTS
def graphs(t_list, function, fourier, K, azero, a, b, funcLabel):
    k_list = np.arange(1,K+1,1)
    
    # FUNCTION
    plt.title("Graph of f(t) and its fourier series")
    plt.plot(t_list, function(t_list), label="f(t) = " + str(funcLabel))
    plt.plot(t_list, fourier, label="Fourier series of f(t) = " + str(funcLabel))
    plt.legend(loc=1)
    plt.xlabel("t")
    plt.ylabel("f(t)")
    plt.savefig("images/exercise1_series_" + str(funcLabel) + ".png", dpi=300)
    plt.show()
    
    # COEFFICIENTS
    plt.title("The coefficients of the fourier series of \nf(t) = " + str(funcLabel))
    plt.scatter(0, azero, color='blue', s=20, label="a$_n$, n $\geq$ 0", marker="x")
    plt.scatter(k_list, a, color='blue', s=20, marker="x")
    plt.scatter(k_list, b, color='red', label="b$_n$, n $\geq$ 1", marker="+")
    plt.legend(loc=1)
    plt.xlim([-.5,K+0.5])
    plt.xlabel("n")
    plt.ylabel("a$_n$, b$_n$")
    plt.savefig("images/exercise1_coefficients_" + str(funcLabel) + ".png", dpi=300)
    plt.show()

def coefficients(K, azero, a, b, funcLabel):
    print("--- Fourier Series coefficients for f(t) = " + str(funcLabel) + "---")
    print("a0 = "+str(azero))
    
    for i in range(1, K+1):
        print("a"+str(i)+" = "+str(a[i-1]))
    print("-")    
    for i in range(1, K+1):
        print("b"+str(i)+" = "+str(b[i-1]))
    
t_list = np.linspace(0,2*np.pi,1000)
K = 10

#################
# f(t) = sin(t) #
#################

# DEFINE FUNCTION
def function(t):
    return np.sin(t)

# GET VALUES
azero, a, b, fourier = fourierSeries(function, 1, 15, K, t_list)

# GRAPHS
graphs(t_list, function, fourier, K, azero, a, b, "sin(t)")

# COEFFICIENTS
coefficients(K, azero, a, b, "sin(t)")

#######################################
# f(t) = cos(t) + 3cos(2t) - 4cos(3t) #
#######################################

# DEFINE FUNCTION
def function(t):
    return np.cos(t) + 3*np.cos(2*t) - 4*np.cos(3*t)

# GET VALUES
azero, a, b, fourier = fourierSeries(function, 1, 15, K, t_list)

# GRAPHS #
graphs(t_list, function, fourier, K, azero, a, b, "cos(t) + 3cos(2t) - 4cos(3t)")

# COEFFICIENTS
coefficients(K, azero, a, b, "cos(t) + 3cos(2t) - 4cos(3t)")

#######################################
# f(t) = sin(t) + 3sin(3t) + 5sin(5t) #
#######################################

# DEFINE FUNCTION
def function(t):
    return np.sin(t) + 3*np.sin(3*t) + 5*np.sin(5*t)

# GET VALUES
azero, a, b, fourier = fourierSeries(function, 1, 20, K, t_list)

# GRAPHS #
graphs(t_list, function, fourier, K, azero, a, b, "sin(t) + 3sin(3t) + 5sin(5t)")

# COEFFICIENTS
coefficients(K, azero, a, b, "sin(t) + 3sin(3t) + 5sin(5t)")

#######################################
# f(t) = sin(t) + 2cos(3t) + 3sin(5t) #
#######################################

# DEFINE FUNCTION
def function(t):
    return np.sin(t) + 2*np.cos(3*t) + 3*np.sin(5*t)

# GET VALUES
azero, a, b, fourier = fourierSeries(function, 1, 20, K, t_list)

# GRAPHS #
graphs(t_list, function, fourier, K, azero, a, b, "sin(t) + 2cos(3t) + 3sin(5t)")

# COEFFICIENTS
coefficients(K, azero, a, b, "sin(t) + 2cos(3t) + 3sin(5t)")
