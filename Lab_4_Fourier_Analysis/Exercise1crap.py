# -*- coding: utf-8 -*-
"""
Created on Tue Mar 28 16:28:17 2023

@author: wattersb
"""

def ak(f, omega, n, K): # FUNCTION, NATURAL FREQUENCY, NUMBER OF STEPS, K STEPS
    ak = np.array([])
    T = 2*np.pi/omega
    
    def akcos(x):
        return f(x) * np.cos(x * k * omega)
    
    for k in range(1,K):
        a = 2/T * Simpson(akcos, 0, T, n)
        ak = np.append(ak, a)
        
    return ak
    
def bk(f, omega, n, K): # FUNCTION, NATURAL FREQUENCY, NUMBER OF STEPS, K STEPS
    bk = np.array([])
    T = 2*np.pi/omega
    
    def bksin(x):
        return f(x) * np.sin(x * k * omega)
    
    for k in range(1,K):
        b = 2/T * Simpson(bksin, 0, T, n)
        bk = np.append(bk, b)
        
    return bk