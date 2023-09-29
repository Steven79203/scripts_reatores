#!/usr/bin/python 
"""
 (a) Para uma operação adiabática e com temperatura inicial de 278 K, 
     plote T, I', rg, −rS, CC e CS em função do tempo até 300 horas. Discuta as tendências observadas.

 (b) Repita (a) e aumente a temperatura inicial com incrementos de 10°C até atingir 330 K e descreva 
     o que você encontrou. Plote a concentração de células após 24 horas de operação em função da temperatura de alimentação.

 (c) Qual a área de troca térmica que deveria ser adicionada para maximizar o número total de 
     células ao final de 24 horas de operação? Para uma temperatura inicial de 310 K e uma 
     temperatura constante de fluido refrigerante de 290 K, qual seria a concentração celular após 24 horas?

"""
# Imports 

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

# == Funções auxiliares == #
# == Velocidade de reação == #
def u_mi(T,Cs):
    c = Cs/(Km+Cs)
    return uMax * I_l(T) * c

def I_l(T):
    a = 0.0038*T*np.exp(21.6 - 6700/T)
    b = 1 + np.exp(153 - 48000/T)
    return a/b

def Qrem(T):
    return U*A*(T-Ta)
    #return mc*Cpc*(T-Ta)*(1-np.exp(1-(U*A)/(mc*Cpc)))

def rC(T,Cs):
    return u_mi(T,Cs) * Cc

def rS(rc):
    return -1/Ycs * rc

def edo(t,F):
    Cs = F[0]
    Cc = F[1]
    T  = F[2]
    dCcdt = u_mi(T,Cs) * Cc
    #dCcdt  = uMax * Cs/(Km+Cs) * Cc
    dCsdt = rS(dCcdt)
    dTdt  = ((dCcdt*V)*(-dHrx) - Qrem(T))/(V*Cps*rho)
    return [dCsdt,dCcdt,dTdt]

# == Dados do problema == #

V    = 25   # dm³
Ycs  = 0.8  # g cel/g subst
Km   = 5.0  # g/dm³
uMax = 0.5  # h^-1
Cps  = 5    # J/g/K
mc   = 100  # kg/h
rho  = 1000 # g/dm³
dHrx = -2e4 # J/g cel
Cpc  = 4.2  # J/g.K
U    = 5e4  # J/(h.K.m²)

A    = 0.31370*1.5 # m²
Ta   = 290   # K

# Condições iniciais do problema
Cc0 = 0.1 # g/dm³
Cs0 = 300 # g/dm³
#T0 = 273+35 # K
T0  = 310 #K

res = solve_ivp(edo,(0,24),(Cs0,Cc0,T0),method='Radau')

t = res.t
Cs = res.y[0]
Cc = res.y[1]
T  = res.y[2]

# Plot
fig1, ax1 = plt.subplots()
ax1.plot(res.t, Cs,label='Cs')
ax1.plot(res.t, Cc,label='Cc')
ax1.set_title('Concentrações (g/dm³)')
ax1.set_ylabel('C (g/dm³)')
ax1.set_xlabel('tempo (h)')
ax1.legend()
ax1.grid()

fig2,ax2 = plt.subplots()
ax2.plot(res.t, T, label='T')
ax2.set_title('Temperatura (K)')
ax2.set_ylabel('T (K)')
ax2.set_xlabel('tempo (h)')
ax2.legend()
ax2.grid()

fig3, ax3 = plt.subplots()
ax3.plot(t, rC(T,Cs),label='rC')
ax3.plot(t, rS(rC(T,Cs)),label='-rS')
ax3.set_title('Taxas de reação (g/dm³h)')
ax3.set_ylabel('-rx(g/dm³.h)')
ax3.grid()
ax3.legend()

fig4, ax4 = plt.subplots()
ax4.plot(t, I_l(T),label="Fator I'")
ax4.set_xlabel('tempo (h)')
ax4.set_ylabel("I'")
ax4.set_title("Adimensional I'")
ax4.grid()
ax4.legend()

fig5, ax5 = plt.subplots()
ax5.plot(res.t, Qrem(T),label='Qrem')
ax5.set_title('Calor removido(J/h)')
ax5.set_ylabel('Q (J/h)')
ax5.set_xlabel('tempo (h)')
ax5.legend()
ax5.grid()

plt.show()
