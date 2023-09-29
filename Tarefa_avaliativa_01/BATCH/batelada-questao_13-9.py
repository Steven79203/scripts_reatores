#!/usr/bin/python

"""
A reação do Problema P11-3A será conduzida em um reator batelada de 10 dm3. 
Plote e analise a temperatura e a concentração de A, B e C em função do tempo 
para os seguintes casos:

(a)Operação adiabática.

(b)Valores de UA de 10.000, 40.000 e 100.000 J/(min ∙ K).

(c)Use UA = 40.000 J/(min ∙ K) e diferentes temperaturas iniciais para o reator.
"""

# A + B -> C
# Batelada

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# Funções auxiliares
def rx(T,X):
    return -(ki(k1,T)*Ca0*(1-X)*Ca0*(theta_b - X))

def ki(k,T):
    return k*np.exp(E/R * (1/Tr - 1/T))

def dHrx(T):
    return deltaHrx + dCp*(T - Tr)

def CiCpi(X):
    Ca = Ca0*(1-X)
    Cb = Ca0*(theta_b - X)
    Cc = Ca0*X
    return Ca*Cpa + Cb*Cpb + Cc*Cpc

def Qgen(T,X):
    return (-rx(T,X)*V)*(-dHrx(T))

def Qrem(T,UA):
    return UA * (T - Ta)

def edo(t,F,UA):
    X = F[0]
    T = F[1]
    dXdt = -rx(T,X)/Ca0
    dTdt = (Qgen(T,X)-Qrem(T,UA))/(V*CiCpi(X))
    return [dXdt, dTdt]

# Dados do problema 

a = 1
b = 1
c = 1

Ca0 = 0.1 # mol/dm³
Cb0 = 0.1 # mol/dm³

theta_b = Cb0/Ca0

V = 10 # dm³

R = 1.987 # cal/mol.K

T0 = 27+273 # K
Tr = 273    # K
Ta = 20+273    # K

Ha = -20000 # cal/mol
Hb = -15000 # cal/mol
Hc = -41000 # cal/mol

Cpa = 15  # cal/mol.K
Cpb = Cpa # cal/mol.K
Cpc = 30  # cal/mol.K


# Fator de conversão J -> cal
f = 0.2390

UA1 = 10000/60  * f # cal/s.K
UA2 = 40000/60  * f # cal/s.K
UA3 = 100000/60 * f # cal/s.K
Ux = [UA1, UA2, UA3]

E = 10000 # cal/mol

k1 = 0.01  # dm³/mol.s at 300K

# Dados adicionais

dCp = (c/a)*Cpc - (b/a)*Cpb - Cpa # cal/mol.K

deltaHrx = Hc - Hb - Ha # cal/mol


Tx = [27+273, 35+273, 45+273]

for Ti in Tx:
    TT = tuple([Ti])
    UU = tuple([UA2])
    res = solve_ivp(edo, (0,120), (0,Ti),method="Radau", args=UU)
    fig1, ax1 = plt.subplots()
    ax1.plot(res.t, Ca0*(1-res.y[0]),label='Ca')
    ax1.plot(res.t, Ca0*(theta_b-res.y[0]),label='Cb')
    ax1.plot(res.t, Ca0*res.y[0],label='Cc')
    ax1.set_xlabel('tempo(s)')
    ax1.set_ylabel('C (mol/dm³)')
    ax1.set_title(f'Concentração - T = {Ti:2.2f} K')
    ax1.legend()
    ax1.grid()
    
    fig2, ax2 = plt.subplots()
    ax2.plot(res.t, res.y[1],label='T',color='r')
    ax2.set_xlabel('tempo(s)')
    ax2.set_ylabel('T (K)')
    ax2.set_title(f'Temperatura - T = {Ti:2.2f} K')
    ax2.legend()
    ax2.grid()

plt.show()
