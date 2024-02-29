# -*- coding: utf-8 -*-
"""
Created on Tue Jan 24 16:36:08 2023

@author: 
"""

import numpy as np
import matplotlib.pyplot as plt

x = np.arange(-5.0, 5.0, 0.2)
a = 1
b = 1
c = -6

plt.plot(x, a * x * x + b * x + c)
plt.plot(x, 0.0 * x)
plt.show()

#fibonacci

def fib(n):
    '''Return a Fibonacci series up to n'''
    a, b = 0, 1
    result=[0]
    while b < n:
        result.append(b)
        a, b = b, a+b
    return result

fib_series = fib(100)
print(fib_series)
