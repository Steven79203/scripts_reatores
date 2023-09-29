#!/usr/bin/python

import numpy as np
from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# Problema adaptado
# Propilenoglicol sintetizado em um PFR com trocador de calor

# == FUNÇÕES AUXILIARES == #
def ki(T):
    return A * np.exp(-E/(R * T))

def rt(T,X):
    return -ki(T)*Ca0*(1-X)

def deltaHrx(T):
    return dHrx + dCp*(T-Tr)

def FiCpi(X):

    Fa = Fa0*(1-X)
    Fb = Fa0*(1-X)
    Fc = Fa0*X
    return Fa*Cpa + Fb*Cpb + Fc*Cpc + Fi0*Cpi

def Qrem(T,Ta):
    return Ua*((T+460.67)-Ta)


# == FUNÇÃO PRINCIPAL == #
def pfr(V,F):
    X  = F[0]
    T  = F[1]
    Ta = F[2]
    dXdV  = -rt(T+460.67,X)/Fa0
    dTdV  = ((-rt(T+460.67,X)*(-deltaHrx(T+460.67))) - Qrem(T,Ta))/(FiCpi(X))
    dTadV = (Ua*((T+460.67)-Ta))/(m*Cpf)

    return [dXdV,dTdV,dTadV]

a = 1
b = 1
c = 1

# Capacidade calorífica das espécies
Cpa = 35    # BTU/lbmol.F
Cpb = 19.5  # BTU/lbmol.F 
Cpi = 18    # BTU/lbmol.F  
Cpc = 46    # BTU/lbmol.F

# Ta trocador de calor
Ua  = 75   # BTU/ft³.h.R
m   = 1102 # lb/h
Cpf = 6.68 # BTU/lb.R

# Entalpias padrão de formação
Ha = -66600  # BTU/lbmol
Hb = -123000 # BTU/lbmol
Hc = -226000 # BTU/lbmol

# Coef de Arrhenius
A = 16.96e12
E = 32400
R = 1.987

# deltaHrx
dHrx = c/a*Hc - b/a*Hb - Ha

# deltaCpi
dCp  = c/a*Cpc - b/a*Cpb - Cpa

# Vazão molar (lbmol/h)
Fa0 = 43.04/10  # Óxido de propileno
Fb0 = 71.87/10  # Metanol
Fi0 = 802.8/10  # Água

# Vazão volumétrica (ft³/h)
va0 = 46.62
vb0 = 233.9
vm0 = 46.62
v0  = (va0 + vb0 + vm0)

# Concentrações das espécies
Ca0 = Fa0/v0  #lbmol/ft²
Cb0 = Fb0/v0  #lbmol/ft²
Ci0 = Fi0/v0  #lbmol/ft²

# Limite de temperatura 584.67 R

# Temperaturas (R)
Ta0 = 545*1  # Temperatura trocador
Tr  = 528       # Temperatura de referência para as entalpia
Ti0 = 75        # Temperatura de alimentação

# Volume de integração
Vf = 100*10 #ft³

res = solve_ivp(pfr,(0,Vf),(0,Ti0,Ta0),method="Radau")

print(res)

V  = res.t
X  = res.y[0]
T  = res.y[1]
Ta = res.y[2]

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()

ax1.plot(V,  T+460.67, label='T')
ax1.plot(V, Ta, label='Ta')
ax1.set_xlabel('V (ft³)')
ax1.set_ylabel('T(R)')
ax1.set_title("Temperaturas")
ax1.legend()
ax1.grid()

ax2.plot(V, X, label='X')
ax2.set_xlabel('V (ft³)')
ax2.set_ylabel('X (lbmol/lbmol)')
ax2.set_title('Conversão')
ax2.legend()
ax2.grid()

plt.show()
