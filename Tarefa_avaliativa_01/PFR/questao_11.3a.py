#!/usr/bin/python 

# Problema 11.3a

# Reação A + B -> C

# Reação adiabática em reator de escoamento contínuo. Uma alimentação equimolar de A e B entra a 27ºC e vazão volumétrica de 2 dm³/s com Ca0 = 0.1 kmol/m³

# a) Grafique e analise a conversão e a temperatura como função do volume para um PFR até X=0.85

# b) Máxima temperatura de entrada até o ponto de ebulição 550K

# Imports 

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# == FUNÇÕES AUXILIARES == #
def dHrx(T):
    return deltaHrx + dCp * (T - Tr)

def ki(k,T):
    return k * np.exp(E/R * (1/Tr - 1/T))

def rx(X,T):
    Ca = Ca0 * (1-X)
    Cb = Ca0 * (theta_b - X)
    return -(ki(k,T) * Ca * Cb)

def Qrem(T):
    return UA * (T - Ta)

def pfr(t,F):
    X = F[0]
    T = F[1]

    dXdV = -rx(X,T)/Fa0
    dTdV = (-rx(X,T)*(-dHrx(T)) - Qrem(T))/(Fa0 * (theta_cpi + X*dCp))

    return [dXdV, dTdV]

# == DADOS DO PROBLEMA == #

a = 1
b = 1
c = 1

Ca0 = 0.1 #mol/dm³
Cb0 = Ca0 #mol/dm³
v0  = 2   #dm³/s
Fa0 = Ca0*v0
Fb0 = Fa0
theta_b = Cb0/Ca0

Ha = -20000 #cal/mol
Hb = -15000 #cal/mol
Hc = -41000 #cal/mol
deltaHrx = (c/a)*Hc - (b/a)*Hb - Ha

Cpa = 15 #cal/mol.K
Cpb = Cpa
Cpc = 30 #cal/mol.K
dCp = Cpc - Cpb - Cpa
theta_cpi = Cpa + theta_b*Cpb

k = 0.01 # dm³/mol.S
E = 10000 # cal/mol
R = 1.987 #cal/mol.K

UA = 1000 #cal/dm³.s.K

T0 = 273+77
Tr = 273 
Ta = 273+30

res = solve_ivp(pfr, (0,100),(0,T0),method="Radau")

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()

ax1.plot(res.t, res.y[0],label='X')
ax1.set_xlabel( 'Volume(dm³)')
ax1.set_ylabel('X (mol/mol)')
ax1.set_title('Conversão')
ax1.grid()
ax1.legend()

ax2.plot(res.t, res.y[1],label='T',color='r')
ax2.set_xlabel('Volume(dm³)')
ax2.set_ylabel('T (K)')
ax2.set_title('Temperatura')
ax2.grid()
ax2.legend()

ax3.plot(res.t, Qrem(res.y[1]),label='Q',color='r')
ax3.set_xlabel('Volume(dm³)')
ax3.set_ylabel('T (K)')
ax3.set_title('Calor (cal/dm³.s)')
ax3.grid()
ax3.legend()

plt.show()
