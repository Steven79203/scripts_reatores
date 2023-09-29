#!/usr/bin/python 

# Lista 01 - batelada - P11.10
# Estequiometria
# 2A -> B + C


# Problemas !
# > Reação muito lenta (dias)
# > Exige altas concentrações de reagentes para ser viável

# > Acarreta em aumento de volume reacional e consequentemente dos custos 
#   iniciais do reator

# > Aumentar a temperatura aumenta a taxa de reação, mas pode ser perigoso
#   dependendo da natureza dos reagentes (temp de ebulição, pressão interna no
#   reator, etc.

# Imports 
import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# == Funções auxiliares == #

def ki(T):
    return k * np.exp(E/R * (1/360 - 1/T))

def rate(T,X):
    return -ki(T) * Ca0**2 * (1-X)**2

def Qgen(T,X):
    return (-rate(T,X))*(-dHrx(T,X))

def dHrx(T,X):
    return deltaHrx + dCp*(T - Tr)

def CiCpi(X):
    Ca = Ca0**2 * (1-X)**2
    Cb = Ca0*(b/a*X)
    Cc = Ca0*(c/a*X)

    return Cpa*Ca + Cpb*Cb + Cpc*Cc + Cpi*Ci0

# == Função principal == #

def edo(t,F):
    X = F[0]
    T = F[1]
    
    dXdt = (-rate(T,X))/Ca0
    dTdt = Qgen(T,X)/CiCpi(X)

    return [dXdt,dTdt]


# Dados
# Concentrações
Ca0 = 4*1000 # mol/dm³
Ci0 = 4*1000 # mol/dm³

# Entalpia de reação
deltaHrx = -10000 # cal/mol

# Temperatura inicial e referência
Ti = 273+70 # K
Tr = 298    # K

# Coeficientes estequimétricos
a,b,c = 2,1,1

# Calor específico (cal/mol.K)
Cpa = 18 
Cpb = 9  
Cpc = 9
Cpi = 15 

dCp = (c/a)*Cpc + (b/a)*Cpb - Cpa

# Velocidade específica de reação
k = 10**-6 #dm³/mol.min

# Energia de ativação
E = 6000  # cal/mol
R = 1.987 # cal/K.mol

# Tempo final (minutos)
tf = 60

sol = solve_ivp(edo, (0,tf),(0,Ti))

# Plot 
fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

ax1.plot(sol.t, sol.y[0], label='Conversão (X)')
ax1.grid()
ax1.legend()

ax2.plot(sol.t, sol.y[1], label='Temperatura(K)')
ax2.grid()
ax2.legend()

plt.show()
