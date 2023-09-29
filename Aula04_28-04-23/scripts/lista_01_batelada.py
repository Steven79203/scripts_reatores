#!/usr/bin/python 

# Bibliotecas importadas 
from scipy.integrate import solve_ivp
import numpy as np
import matplotlib.pyplot as plt


# == Funções auxiliares == #
# Velocidade específica da reação
def kin(k,E,T):
    return k * np.exp(E/R * (1/Ti - 1/T))

# Entalpia de geração
def dHrx(T):
    return deltaHrx + dCp*(T - Tr)

# Taxa de consumo
def rate(X,T):
    Ca = Ca0 * (1-X)
    Cb = Ca0 * (theta_b - X)
    Cc = Ca0 * X
    k_1 = kin(k1,E1,T)
    k_2 = kin(k2,E2,T)
    return -(k_1 * Ca**0.5 * Cb**0.5 - k_2*Cc)

# Calor removido
def Qrem(T):
    return UA * (T - Ta)

# Calor gerado
def Qgen(X,T):
    return  (-rate(X,T)*V*-dHrx(T))

# Função EDO
def edo(t,F):
    dXdt = (-rate(F[0],F[1]))/Ca0
    dTdt = ((Qgen(F[0],F[1]) - Qrem(F[1]))/NiCpi)
    return [dXdt,dTdt]


# == Dados do problema == #
# Velocidades específicas de reação (373K)
k1 = 2e-3 # s^-1
k2 = 3e-5 # s^-1

# Energia padrão de formação 
E1 = 100*1000 # J/mol
E2 = 150*1000 # J/mol

# Volume (dm³)
V = 10

# Concentração inicial dos reagentes
Ca0 = 0.1  # mol/dm³
Cb0 = 0.125 # mol/dm³

# Mols
Na0 = Ca0*V
Nb0 = Cb0*V

# Capacidades caloríficas dos componentes
Cpa = 25 # J/mol.K
Cpb = 40 # J/mol.K
Cpc = 40 # J/mol.K

# NiCpi
NiCpi = Na0*Cpa + Nb0*Cpb

# Entalpia padrão de formação 
deltaHrx = -40000.00 # J/mol.A

# Constante dos gases
R = 8.3144 # J/mol.K

# Temperatura inicial
Ti = 373   # K

# Coeficiente acoplado global de transferêcia de calor
UA = 2424 # J/s.ºC

# Temperatura de referência (trocador de calor)
Ta = 370

# Variação da capacidade calorífica
dCp = Cpc - Cpb - Cpa

# Temperatura de referência
Tr = 293

# == Dados adicionais == #
# Theta (Nx0/Na0)
theta_a = 1
theta_b = Cb0/Ca0
theta_c = 0

# Intervalo de integração
tf = 60*60

# Solver 
sol = solve_ivp(edo,(0,tf), (0,Ti),method="Radau")

# Plot dos dados
fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()

ax1.plot(sol.t, sol.y[1],label='Temperatura (K)')
ax1.set_title('Temperatura')
ax1.set_xlabel('t(seg)')
ax1.set_ylabel('T (K)')

ax2.plot(sol.t, sol.y[0])
ax2.set_title('Conversão')
ax2.set_xlabel('t(seg)')
ax2.set_ylabel('X')

# Concentração das espécies
ax3.plot(sol.t, Ca0*(1-sol.y[0]), label='Ca')
ax3.plot(sol.t, Ca0*(theta_b - sol.y[0]), label='Cb')
ax3.plot(sol.t, Ca0*sol.y[0], label='Cc')
ax3.set_title('Concentração das espécies')
ax3.set_xlabel('tempo(s)')
ax3.set_ylabel('C (mol/dm³)')

ax1.legend()
ax2.legend()
ax3.legend()

ax1.grid()
ax2.grid()
ax3.grid()

plt.show()
