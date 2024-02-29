# -*- coding: utf-8 -*-
"""
Created on Tue Apr  4 14:13:48 2023

@author: 
"""
import numpy as np
import matplotlib.pyplot as plt

# FUNDANMENTAL FREQUENCY
omega = 1

# TAU FOR RECTANGULAR WAVE
tau = 1

# SIMPSONS RULE FOR INTEGRATION
def Simpson(f, a, b, n): # FUNCTION, LOWER LIMIT, UPPER LIMIT, NUMBER OF STEPS
    
    h = (b-a) / n
    
    # THE SUMMATION TERMS
    sigma1 = 0
    sigma2 = 0
    
    j = 1
    while j <= n/2 - 1:
        x2j = a + 2*j*h
        sigma1 = sigma1 + f(x2j)

        j += 1

    j = 1
    while j <= n/2:
        xj = a + (2*j-1)*h
        sigma2 = sigma2 + f(xj)
        j += 1
    
    return h/3 * (f(a) + 2*sigma1 + 4*sigma2 + f(b))


#############################################
#                                           #
#   DEFINE FUNCTIONS TO FIND COEFFICIENTS   #
#                                           #
#############################################


def a0(f, omega, n): # NATURAL FREQUENCY, FUNCTION, NUMBER OF STEPS
    T = 2*np.pi/omega
    return 1/T * Simpson(f, 0, T, n)

def ak(f, omega, n, k): # NATURAL FREQUENCY, FUNCTION, NUMBER OF STEPS, k STEPS
    T = 2*np.pi/omega
    
    # EMPTY ARRAY OF ak VALUES
    ak = np.array([])
    
    # THE FUNCTION TIMES THE COSINE TERM
    def fcos(x):
        return f(x)*np.cos(x*omega*i)
    
    # CALCULATE THE ak TERMS
    for i in range(1, k+1):
        a = Simpson(fcos, 0, T, n)
        ak = np.append(ak, a)
            
    ak = 2/T * ak
    return ak

def bk(f, omega, n, k): # NATURAL FREQUENCY, FUNCTION, NUMBER OF STEPS, k STEPS
    T = 2*np.pi/omega
    
    # EMPTY ARRAY OF bk VALUES
    bk = np.array([])
    
    # THE FUNCTION TIMES THE SINE TERM
    def fsin(x):
        return f(x)*np.sin(x*omega*i)
    
    # CALCULATE THE bk TERMS
    for i in range(1, k+1):
        b = Simpson(fsin, 0, T, n)
        bk = np.append(bk, b)

    bk = 2/T * bk
    return bk

def fourierSeries(f, omega, n, k, t): # NATURAL FREQUENCY, FUNCTION, NUMBER OF STEPS, k STEPS, t VALUE
    
    # CALCULATE THE COEFFICIENTS
    azero = a0(f, omega, n)
    a = ak(f, omega, n, k)
    b = bk(f, omega, n, k)

    # CALCULATE THE SERIES
    fourier = 0
    for n in range(1, k+1):
        fourier += a[n-1]*np.cos(n*t) + b[n-1]*np.sin(n*t)
    fourier += azero
    
    return azero, a, b, fourier


###################
#                 #
#   SQUARE WAVE   #
#                 #
###################


def square(t_list):
    
    step_list = np.array([])
    
    # USED TO ENSURE THE FUNCTION IS PLOTTED CORRECTLY OVER ALL RANGES
    def squareFunc(t_list):
        theta = t_list * omega
        
        # KEEP THETA WITHIN [0, 2pi]
        while np.abs(theta) > 2*np.pi:
             theta = theta - 2*np.pi * np.abs(theta)/theta
        
        # ASSIGN VALUES
        if (theta >= 0 and theta <= np.pi) or (theta < -np.pi and theta >= -2*np.pi):
            return 1
        else:
            return -1 
    
    # WORKS SLIGHTLY DIFFERENTLY FOR NUMBERS AND ARRAYS - PYTHON WON'T "LOOP" OVER NUMBERS
    if type(t_list) == float or type(t_list) == int:
        step_list = np.append(step_list, squareFunc(t_list))
    else:
        for t in t_list:
            step_list = np.append(step_list, squareFunc(t))
    
    return step_list


# GET VALUES
t_list = np.linspace(0,2*np.pi,2000)
k_values = np.array([1,2,3,5,10,20,30])
azero, a, b, fourier = fourierSeries(square, 1, 1000, k_values[-1], t_list)

##########
# GRAPHS #
##########

# TOGETHER
plt.title("Graph of the square wave with its Fourier series for\nvarious values of n superimposed")
plt.plot(t_list, square(t_list), label="Square function")
for i in k_values:
    plt.plot(t_list, fourierSeries(square, 1, 1000, i, t_list)[3], label="n=" + str(i))

plt.legend(loc=1)
plt.xlabel("t")
plt.ylabel("f(t)")
plt.savefig("images/exercise2_square.png", dpi=300)
plt.show()

# DIFFERENT K VALUES
for i in k_values:
    plt.title("Graph of the square wave with its Fourier series for n = "+str(i))
    plt.plot(t_list, square(t_list), label="Square function")
    plt.plot(t_list, fourierSeries(square, 1, 1000, i, t_list)[3], label="Fourier series with n=" + str(i)  )
    plt.legend(loc=1)
    plt.xlabel("t")
    plt.ylabel("f(t)")
    plt.savefig("images/exercise2_square_k" + str(i) + ".png", dpi=300)
    plt.show()
    
# OUTSIDE OF [0, 2pi], k = 30
t_list = np.linspace(-10,10,5000)

plt.title("Graph of the square wave with its Fourier series for n = 30")
plt.plot(t_list, square(t_list), label="Square function")
plt.plot(t_list, fourierSeries(square, 1, 1000, 30, t_list)[3], label="Fourier series with n=30")
plt.legend(loc=1)
plt.xlabel("t")
plt.ylabel("f(t)")
plt.savefig("images/exercise2_square_long.png", dpi=300)
plt.show()

# PRINT COEFFICIENTS
print("--- Fourier Series coefficients for the square wave function ---")
print("a0 = "+str(azero))

for i in range(1, k_values[-1]+1):
    print("a"+str(i)+" = "+str(a[i-1]))
print("-")    
for i in range(1, k_values[-1]+1):
    print("b"+str(i)+" = "+str(b[i-1]))


########################
#                      #
#   RECTANGULAR WAVE   #
#                      #
########################


def rectangular(t_list):
    
    # USED TO ENSURE THE FUNCTION IS PLOTTED CORRECTLY OVER ALL RANGES
    def rectangularFunc(t_list):
        theta = t_list * omega
        
        # KEEP THETA WITHIN [0, 2pi]
        while np.abs(theta) > 2*np.pi:
             theta = theta - 2*np.pi * np.abs(theta)/theta
        
        # ASSIGN VALUES
        if (theta >= 0 and theta <= omega * tau) or (theta >= -2*np.pi and theta <= -2*np.pi + omega * tau):
            return 1
        else:
            return -1 
        
    step_list = np.array([])
    
    # WORKS SLIGHTLY DIFFERENTLY FOR NUMBERS AND ARRAYS - PYTHON WON'T "LOOP" OVER NUMBERS
    if type(t_list) == float or type(t_list) == int:
        step_list = np.append(step_list, rectangularFunc(t_list))
    else:
        for t in t_list:
            while np.abs(t) > 2* np.pi:
                t = t - 2*np.pi * np.abs(t)/t
            step_list = np.append(step_list, rectangularFunc(t))
    
    return step_list


# GET VALUES
t_list = np.linspace(0,2*np.pi,2000)
k_values = np.array([1,2,3,5,10,20,30])
azero, a, b, fourier = fourierSeries(rectangular, 1, 1000, k_values[-1], t_list)

##########
# GRAPHS #
##########

# TOGETHER
plt.title("Graph of the rectangular wave with its Fourier series for\nvarious values of n superimposed")
plt.plot(t_list, rectangular(t_list), label="Rectangular function")
for i in k_values:
    plt.plot(t_list, fourierSeries(rectangular, 1, 1000, i, t_list)[3], label="n=" + str(i))
    
plt.legend(loc=1)
plt.xlabel("t")
plt.ylabel("f(t)")
plt.savefig("images/exercise2_rect.png", dpi=300)
plt.show()

# DIFFERENT K VALUES
for i in k_values:
    plt.title("Graph of the rectangular wave with its Fourier series for n = "+str(i))
    plt.plot(t_list, rectangular(t_list), label="Rectangular function")
    plt.plot(t_list, fourierSeries(rectangular, 1, 1000, i, t_list)[3], label="Fourier series with n=" + str(i))
    plt.legend(loc=1)
    plt.xlabel("t")
    plt.ylabel("f(t)")
    plt.savefig("images/exercise2_rect_k" + str(i) + ".png", dpi=300)
    plt.show()

# OUTSIDE OF [0, 2pi], k = 30
t_list = np.linspace(-10,10,5000)

plt.title("Graph of the rectangular wave with its Fourier series for n = 30")
plt.plot(t_list, rectangular(t_list), label="Rectangular function")
plt.plot(t_list, fourierSeries(rectangular, 1, 1000, 30, t_list)[3], label="Fourier series with n=30")
plt.legend(loc=1)
plt.xlabel("t")
plt.ylabel("f(t)")
plt.savefig("images/exercise2_rect_long.png", dpi=300)
plt.show()

# PRINT COEFFICIENTS
print("--- Fourier Series coefficients for the rectangular wave function ---")
print("a0 = "+str(azero))

for i in range(1, k_values[-1]+1):
    print("a"+str(i)+" = "+str(a[i-1]))
print("-")    
for i in range(1, k_values[-1]+1):
    print("b"+str(i)+" = "+str(b[i-1]))










