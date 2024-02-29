# -*- coding: utf-8 -*-
"""
Created on Wed Apr  5 11:39:47 2023

@author: 
"""
import numpy as np
import matplotlib.pyplot as plt


########################
#                      #
#   DEFINE FUNCTIONS   #
#                      #
########################


#####################
# FOURIER TRANSFORM #
#####################

# REAL PART
def fReal(f, N, n, h):
    
    sigma = 0
    for m in range(N):
        sigma += f(m*h) * np.cos(2*np.pi*m*n/N)
    
    return sigma
    
# IMAGINARY PART
def fImag(f, N, n, h):
    
    sigma = 0
    for m in range(N):
        sigma += -f(m*h) * np.sin(2*np.pi*m*n/N)
    
    return sigma


##################
# BACK TRANSFORM #
##################

def fmReal(Freal, Fimag, N, m):
    sigma = 0
    
    for n in range(N):
        sigma += Freal(n) * np.cos(2*np.pi*m*n/N) - Fimag(n) * np.sin(2*np.pi*m*n/N)
    
    return 1/N * sigma

def fmImag(Freal, Fimag, N, m):
    sigma = 0
    
    for n in range(N):
        sigma += Freal(n) * np.sin(2*np.pi*m*n/N) + Fimag(n) * np.cos(2*np.pi*m*n/N)
    
    return 1/N * sigma


####################
#                  #
#   sin(0.45pit)   #
#                  #
####################


# FUNCTION BEING ANALYSED
def func(t):
    return np.sin(0.45*np.pi*t)

#############
# TRANSFORM #
#############

N = 128 # number of samples
h_list = [0.1, 5/144] # sampling interval
h_labels = ["0.1", "5/144"]

x = 0
while x < len(h_list):
    
    h = h_list[x]
    
    # SAMPLING RATE
    rate = 1/h
    print("The sampling rate is " + str(rate) + " samples per second")
    
    # FUNDAMENTAL FREQUENCY
    omega1 = 2*np.pi/(h*N)
    print("The fundamental frequency is " + str(omega1) + " /s")
    
    # PLOT THE FUNCTION
    t_list = np.linspace(0, (N-1) * h, 1000)
    plt.plot(t_list, func(t_list))
    
    # SAMPLE TIMES
    m = np.arange(0, N, 1)
    tm = m*h
    
    # ADD SAMPLED POINTS
    plt.scatter(tm, func(tm), marker='x', color='black')
    plt.xlabel("t")
    plt.ylabel("f(t)")
    plt.title("The graph of the function sin(0.45$\pi$t) for h = " + h_labels[x])
    plt.savefig("images/exercise3_0.45_graph_"+str(h)+".png", dpi=300)
    plt.show()
    
    real = np.array([])
    imag = np.array([])
    
    n = np.arange(1, N+1, 1)
    for i in n:
        real = np.append(real, fReal(func, N, i, h))
        imag = np.append(imag, fImag(func, N, i, h))
    
    plt.title("The components of the Fourier Transform of sin(0.45$\pi$t) for h = " + h_labels[x])
    plt.plot(n, real, label="Real component of Fourier Transform")
    plt.plot(n, imag, label="Imaginary component of Fourier Transform")
    plt.xlabel("n")
    plt.ylabel("F$_n$")
    plt.legend()
    plt.savefig("images/exercise3_0.45_components_"+str(h)+".png", dpi=300)
    plt.show()
    
    ##################
    # BACK TRANSFORM #
    ##################
    
    def realF(n):
        return fReal(func, N, n, h)
    
    def imagF(n):
        return fImag(func, N, n, h)

    back_real = np.array([])
    back_imag = np.array([])
    
    for i in m:
        back_real = np.append(back_real, fmReal(realF, imagF, N, i))
        back_imag = np.append(back_imag, fmImag(realF, imagF, N, i))
    
    plt.title("The reconstructed function and from the Fourier back-transform  \n for h = " + h_labels[x] + " and the original function")
    plt.xlabel("t")
    plt.ylabel("f(t)")
    plt.plot(t_list, func(t_list), label="Original function", color="blue")
    plt.plot(tm, back_real, label="Fourier back-transform (real)", color="red")
    plt.plot(tm, back_imag, label="Fourier back-transform (imaginary)", color="yellow")
    plt.legend()
    plt.savefig("images/exercise3_0.45_back_"+str(h)+".png", dpi=300)
    plt.show()
    
    plt.title("The difference between the reconstructed and \n original functions for h = " + h_labels[x])
    plt.plot(tm, func(tm) - back_real)
    plt.xlabel("t")
    plt.ylabel("error")
    plt.savefig("images/exercise3_0.45_error_"+str(h)+".png", dpi=300)
    plt.show()
    
    x += 1
    
#################
#               #
#   cos(6pit)   #
#               #
#################


# FUNCTION BEING ANALYSED
def func(t):
    return np.cos(6*np.pi*t)


# LOOP OVER h VALUES
h_list = np.array([0.6, 0.2, 0.1, 0.04])

# INITIAL VARIABLES
N = 32
m = np.arange(0, N, 1)
n = np.arange(1, N+1, 1)

for h in h_list:
    #############
    # TRANSFORM #
    #############
    
    N = 32 # number of samples
    
    # SAMPLING RATE
    rate = 1/h
    print("The sampling rate is " + str(rate) + " samples per second")
    
    # FUNDAMENTAL FREQUENCY
    omega1 = 2*np.pi/(h*N)
    print("The fundamental frequency is " + str(omega1) + " /s")
    
    # PLOT THE FUNCTION
    t_list = np.linspace(0, (N-1) * h, 1000)
    plt.plot(t_list, func(t_list))
    
    # SAMPLE TIMES
    tm = m*h
    
    # ADD SAMPLED POINTS
    plt.scatter(tm, func(tm), marker='x', color='black')
    plt.xlabel("t")
    plt.ylabel("f(t)")
    plt.title("The graph of the function cos(6$\pi$t) for h = " + str(h))
    plt.savefig("images/exercise3_0.6_graph_"+str(h)+".png", dpi=300)
    plt.show()
    
    real = np.array([])
    imag = np.array([])
    
    for i in n:
        real = np.append(real, fReal(func, N, i, h))
        imag = np.append(imag, fImag(func, N, i, h))
    
    plt.title("The components of the Fourier Transform of cos(6$\pi$t) for h = " + str(h))
    plt.plot(n, real, label="Real component of Fourier Transform")
    plt.plot(n, imag, label="Imaginary component of Fourier Transform")
    plt.xlabel("n")
    plt.legend()
    plt.savefig("images/exercise3_0.6_components_"+str(h)+".png", dpi=300)
    plt.show()
    
    ##################
    # POWER SPECTRUM #
    ##################
    
    def powerSpec(func, N, n, h):
        return fReal(func, N, n, h)**2 + fImag(func, N, n, h)**2
    
    plt.title("Power spectrum for h = " + str(h))
    plt.plot(n, powerSpec(func, N, n, h))
    plt.xlabel("n")
    plt.savefig("images/exercise3_0.6_power_"+str(h)+".png", dpi=300)
    plt.show()
    
    ##################
    # BACK TRANSFORM #
    ##################
    
    def realF(n):
        return fReal(func, N, n, h)
    
    def imagF(n):
        return fImag(func, N, n, h)
    
    back_real = np.array([])
    back_imag = np.array([])
    
    for i in m:
        back_real = np.append(back_real, fmReal(realF, imagF, N, i))
        back_imag = np.append(back_imag, fmImag(realF, imagF, N, i))
    
    plt.title("The reconstructed function and from the Fourier back-transform  \n for h = " + str(h) + " and the original function")
    plt.xlabel("t")
    plt.ylabel("f(t)")
    plt.plot(t_list, func(t_list), label="Original function", color="blue")
    plt.plot(tm, back_real, label="Fourier back-transform (real)", color="red")
    plt.plot(tm, back_imag, label="Fourier back-transform (imaginary)", color="yellow")
    plt.legend()
    plt.savefig("images/exercise3_0.6_back_"+str(h)+".png", dpi=300)
    plt.show()
    
    plt.title("The difference between the reconstructed and \n original functions for h = " + str(h))
    plt.plot(tm, func(tm) - back_real)
    plt.xlabel("t")
    plt.ylabel("error")
    plt.savefig("images/exercise3_0.6_error_"+str(h)+".png", dpi=300)
    plt.show()
