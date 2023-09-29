#!/usr/bin/python 

from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt

# Produção de nitroanilina - Acidente de Monsanto
# Sem excesso de reagente 

# Legenda
# A = ONBC
# B = NH3
# C = H20
# D = Nitroanilina

# == Funções Auxiliares == #
# Constante de velocidade 
def kin(T):
    return 0.00017 * np.exp(11273/1.987 * (1/461 - 1/T))

# Taxa de reação
def rate(k,X):
    return -k * C[0]**2 * (1-X)*(theta[1] - 2*X)

# Calor gerado 
def Qgen(T,X):
    return (rate(kin(T),X)*V)*(dHrx)

# Calor removido
def Qrem(T):
    return UA * (T - Ta)

# Sistema de EDOs
def edo(t,F):
    Qgen_i = Qgen(F[1],F[0])
    Qrem_i = 0

    if t < 45:
        Qrem_i = Qgen_i
    elif t >= 45 and t < 55:
        Qrem_i = 0 
    else:
        Qrem_i = Qrem(F[1])

    dTdt = (Qgen_i - Qrem_i)/NiCpi
    dXdt = -rate(kin(F[1]),F[0]) * V / Nx[0]
    
    return (dXdt, dTdt)

# Alimentação inicial normal (kmol)
Nx = np.array([3.17,43,103.6])

# Alimentação extra
#Nx = np.array([9.044, 33.0, 103.7])

# Theta
theta = np.array(Nx/Nx[0])

# Capacidades caloríficas (cal/kmol)
Cpx = np.array([40, 8.38, 18])

# Entalpia padrão de reação (kcal/kmol)
dHrx = -5.9e5

# NiCpi (kcal/K)
NiCpi = np.dot(Nx,Cpx)

# Volume reacional (m³)
V = 5.119
#V = 3.26

# Temp de referência (K)
Ta = 298

# UA (kcal/min Cº)
UA = 35.85 

# Concentrações
C = np.array(Nx/V)

# Solver
sol = solve_ivp(edo,(0,120),(0,448),max_step=0.5)

# Plot
fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
ax1.plot(sol.t,sol.y[0],label='Conversão',color='g')
ax1.set_xlabel('tempo(min)')
ax1.set_ylabel('Conversão (kmol/kmol)')
ax2.plot(sol.t,sol.y[1],label='T (K)',color='b')
ax2.set_ylabel('T (K)')

ax1.legend()
ax1.grid()
ax2.legend()
ax2.grid()

fig3, ax3 = plt.subplots()
ax3.plot(sol.t, Qgen(sol.y[1],sol.y[0]), label='Calor gerado', color='r')
ax3.plot(sol.t, Qrem(sol.y[1]), label='Calor removido', color='m')
ax3.set_xlabel('tempo (min)')
ax3.set_ylabel('Calor (kcal/min)')
ax3.grid()
ax3.legend()

plt.show()
#for t, X, T in zip(sol.t, sol.y[0], sol.y[1]):
#    print(f"{t:^12.5}   |   {X:^12.5}   |   {T:^12.5}")

