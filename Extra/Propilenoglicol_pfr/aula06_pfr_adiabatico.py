#!/usr/bin/python

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def ki(T):
    return A * np.exp(-E/(R * T))

def rt(T,X):
    return -ki(T)*Ca0*(1-X)

def deltaHrx(T):
    return dHrx + dCp*(T-Tr)

def pfr(V,F):
    X  = F[0]
    T  = F[1]
    Ta = F[2]
    dXdV  = -rt(T,X)/Fa0
    dTdV  = ((-rt(T,X)*(-deltaHrx(T))) - Ua*(T-Ta))/(Cp0*Fa0)
    dTadV = (Ua*(Ta - T))/(m*Cpc)

    return [dXdV,dTdV,dTadV]

a = 1
b = 1
c = 1

# Capacidade calorífica das espécies
# BTU/lbmol ºF
Cpa = 35   
Cpb = 18
Cpn = 19.5
Cp0 = 37  

Ua  = 1000 # BTU/ft³.h.R
m   = 1102 # lb/h
Cpc = 6.73 # BTU/lbmol.R

# Entalpias padrão de formação
Ha = -66600   # BTU/
Hb = -123000
Hc = -226000

# Coef de Arrhenius
A = 16.96e12
E = 32400
R = 1.987

# deltaHrx
dHrx = c/a*Hc - b/a*Hb - Ha

# deltaCpi
dCp  = c/a*Cpc - b/a*Cpb - Cpa

# Vazão molar
Fa0 = 43.04*1/100  # Óxido de propileno
Fb0 = 802.8       # Água
Fm0 = 71.87       # Metanol

# Vazão volumétrica
va0 = 46.62
vb0 = 233.9
vm0 = 46.62
v0  = (va0 + vb0 + vm0)

# Concentrações das espécies
Ca0 = Fa0/v0
Cb0 = Fb0/v0
Cm0 = Fm0/v0

# Limite de temperatura 584.67 R

Ta0 = 545  # Temperatura trocador
Tr  = 528  # Temperatura de referência para as entalpia
Ti0 = 535  # Temperatura de alimentação

Vf = 25

res = solve_ivp(pfr,(0,Vf),(0,Ti0,Ta0),method="Radau")

print(res)

V  = res.t
X  = res.y[0]
T  = res.y[1]
Ta = res.y[2]

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

ax1.plot(V,  T, label='T')
ax1.plot(V, Ta, label='Ta')

ax2.plot(V, X, label='X')

plt.legend()
plt.grid()
plt.show()
