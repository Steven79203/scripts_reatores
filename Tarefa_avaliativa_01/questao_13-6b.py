#!/usr/bin/python

# Reação A + B -> C

# -ra = k1.Ca^(1/2).Cb^(1/2) - k2Cc

import numpy as np
from scipy.integrate import solve_ivp
from matplotlib.pyplot as plt

def dHrx(T):
    return deltaHrx + dCp * (T - Tr)

def ki(k,E,T):
    return k * np.exp(E/R*(1/Tr - 1/T))

def rt(X,T):
    return -(ki(k1,E1,T)*Ca0*(1-X)*(theta_b-X)-ki(k2,E2,T)*Ca0*X)

def Qgen(X,T):
    return (-rt(X,T))*V*(-dHrx(T)))

def NiCpi(X):
    Na = Na0/V * (1-X)
    Nb = Na0/V * (theta_b -X)
    Nc = Na0/V * X
    return Na*Cpa + Nb*Cpb + Nc*Cpc

def ode(t,F):
    X = F[0]
    T = F[1]
    dXdt = -rt(X,T)/Ca0
    dTdt = Qgen(X,T)/NiCpi(X)
    return [dXdt, dTdt]

a =  1
b =  1
c =  1

V  = 10   # dm³

k1 = 2e-3 # s^-1
k2 = 3e-5 # s^-1

E1 = 100 # kJ/mol
E1 = 150 # kJ/mol

Ca0 = 0.1   # mol/dm³
Cb0 = 0.125 # mol/dm³

theta_a = 1
theta_b = Cb0/Ca0

deltaHrx = -40000 # J/mol.A

Cpa = 25 # J/mol.K
Cbp = 25 # J/mol.K
Cpc = 40 # J/mol.K

dCp = (c/a)*Cpc - (b/a)*Cpb - (a/a) 

Tr = 298   # K

R  = 8.314 # J/K.mol


