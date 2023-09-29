#!/usr/bin/python 

# Lista 01 - Exercício 04 - P12-6

# Estequiometria
# a + b -> 2c

import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp

# == Equações auxiliares == #
# Velocidade específica
def ki(k,T):
    return k * np.exp(E/R * (1/(80+459.69) - 1/(T+459.69)))

# Taxa de reação
def rt(T,X):
    return -ki(k1,T) * Ca0**2 * (1-X) * (theta_b - X)

# Calor de reação
def dHrx(T):
    return deltaHrx + dCp*(T - Ti)

# Somatório Ci e Cpi
def CiCpi(X):
    Ca = Ca0*(1-X)
    Cb = Ca0*(theta_b - X)
    Cc = Ca0*(2*X)
    return Cpa*Ca + Cpb*Cb + Cpc*Cc

# Calor de reação
def Qr(T,X):
    return (-rt(T,X)*V)*(-dHrx(T))

# Trocador de calor
def Qex(T):
    return U * A * (Tvap - T)

# Calor do agitador
def Wsp():
    return Ws * 2.545 # btu/h

# == EDO == #
def edo(t,F):
    X = F[0]
    T = F[1]
    dXdt = (-rt(T,X)*V)/Na0
    dTdt = (Qex(T) + Wsp() + Qr(T,X))/(V*CiCpi(X))
    return [dXdt, dTdt]

# == Dados == #

V    = 125 # gal 
A    = 10  # ft²
Pvap = 150 # psi
Tvap = 365 # ºF
U    = 150 # btu/h.ft².ºF
Ws   = 25*2.545  # btu/h

# Mols (lbmol)
Na0 = 1000
Nb0 = 1000

# Concentração
Ca0 = Na0/V
Cb0 = Nb0/V

# Theta
theta_a = 1
theta_b = Nb0/Na0
theta_c = 0

# Massa molar (lb/lbmol)
Mma = 128 
Mmb = 94
Mmc = 111

# Massa específica (lbmol/ft³)
rho_a = 63.0
rho_b = 67.2
rho_c = 65.0

# Calor específico (but/lbmol ºF)
Cpa = 51   
Cpb = 44   
Cpc = 47.5 

# DeltaHrx
deltaHrx = 20000 # btu/lbmol

# dCp
dCp = 2*Cpc - Cpb - Cpa

# == Dados cinéticos == #
k1 = 5      # h^-1 (80ºF)
Ti = 80     # ºF
E  = 10000  # but/lbmol
R  = 1.987  # but/(lbmol.R)

tf = 1 # horas

sol = solve_ivp(edo,(0,tf),(0,Ti))

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

ax1.plot(sol.t, sol.y[0], label='Conversão(X)')
ax1.set_xlabel('tempo (h)')
ax1.set_ylabel('X (lbmol/lbmol)')

ax2.plot(sol.t, sol.y[1], label='Temperatura(ºF)')
ax2.set_xlabel('tempo(h)')
ax2.set_ylabel('T (ºF)')

ax1.grid()
ax2.grid()

ax1.legend()
ax2.legend()

plt.show()
