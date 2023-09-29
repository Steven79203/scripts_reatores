#!/usr/bin/python

import numpy as np
from scipy.optimize import fsolve
import matplotlib.pyplot as plt

def ki(T):
    return A * np.exp(-E/(R * T))

def rt(T,X):
    return -ki(T)*(Fa0/v0)*(1-X)

def deltaHrx(T):
    return dHrx + dCp*(T-Tr)

def cstr(x):
    X = x[0]
    T = x[1]
    F = V*(-rt(T,X)) - Fa0*X
    G = - theta_cpi*(T-To) - deltaHrx(T)*X
    return [F,G]

a = 1
b = 1
c = 1

Cpa = 35
Cpb = 18
Cpc = 46
Cpn = 19.5

Ha = -66600
Hb = -123000
Hc = -226000

V   = 40

A = 16.96e12
E = 32400
R = 1.987

UA = 100*40

dHrx = c/a*Hc - b/a*Hb - Ha
dCp  = c/a*Cpc - b/a*Cpb - Cpa

Fa0 = 43.04
Fb0 = 802.8
Fn0 = 71.87

theta_b = Fb0/Fa0
theta_n = Fn0/Fa0
theta_cpi = Cpa + theta_b*Cpb + theta_n*Cpn

va0 = 46.62
vb0 = 233.9
vn0 = 46.62
v0  = va0 + vb0 + vn0

Ta = 545
Tr = 528
To = 535

X0 = 0.9
T0 = 600

res = fsolve(cstr, [X0,T0])

print(res)

Xi = np.linspace(0,1,50)
Ti = np.linspace(300,650,50)

plt.plot(Xi,cstr([Xi,Ti])[0])
plt.show()
