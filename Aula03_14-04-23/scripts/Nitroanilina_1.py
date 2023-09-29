#!/usr/bin/python 

from scipy.integrate import solve_ivp
from sys import argv
import numpy as np
import matplotlib.pyplot as plt


# Produção de nitroanilina - Acidente de Monsanto

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
    return -k * C[0]**2 * (1-X) *(theta[1] - 2*X)

# Calor gerado 
def Qgen(T,X):
    return (rate(kin(T),X)*V)*(dHrx)

# Calor removido
def Qrem(T):
    return UA * (T - Ta)

def Qrem_2(T):
    return mc * Cpfluido * (T-Ta) * (1-np.exp(-UA/mc/Cpfluido))

# == CALLBACKS == #
# Refrigeração falha e não há disco de ruptura (RFND)
def refr_fail_nodisk(t,T,X):
    if t < 45:
        return Qgen(T,X)
    if t >= 45 and t < 55:
        return 0
    if t >= 55:
        return Qrem(T)

# Refrigeração falha, mas há disco de ruptura (RFWD)
def refr_fail_disk(t,T,X):
    if t < 45:
        return Qgen(T,X)
    if t >= 45 and t < 55:
        return 0
    if t >= 55 and t < 114:
        return Qrem(T)
    if t >= 114 and t < 114+2.24:
        return Qrem(T) + (830*538 if disk_ruptured else 0)
    if t >= 114+2.24:
        return Qrem(T)

# Refrigeração nunca falha (RNF)
def refr_nofail_nodisk(T):
    return Qrem(T)

def refr_fail_more(t,T,X):
    if t < 45:
        return Qgen(T,X)
    if t >= 45:
        return (830*538 if disk_ruptured else 0)

callbacks = {
        "rfnd" : refr_fail_nodisk,
        "rfyd" : refr_fail_disk,
        "rnf"  : refr_nofail_nodisk,
        "rfm"  : refr_fail_more
        }

# Sistema de EDOs
def edo(t,F):

    global disk_ruptured
    if F[1] > 538:
        disk_ruptured = True

    Qgen_i = Qgen(F[1],F[0])
    Qrem_i = refr_case(t,F[1],F[0])

    dTdt = (Qgen_i - Qrem_i)/NiCpi
    dXdt = -rate(kin(F[1]),F[0]) * V / Nx[0]
    
    return (dXdt, dTdt)


# Estado do disco
global disk_ruptured, disk_t, disk_i
disk_ruptured = False # Disco rompido

# Alimentação inicial normal (kmol)
Nx0 = np.array([3.17,43,103.6])

# Alimentação extra
Nx = np.array([9.044, 33.0, 103.7])

# Alimentação normal e extra
Nv = np.array([Nx0,Nx])

# Volume reacional (m³)
Vexc  = 5.119
Vnorm = 3.26

# Definir excesso de reagentes e volume
# PADRÃO: Sem excesso e volume normal
try:
    if (int(argv[2]) == 1):
        mols = Nv[1]
        V = Vexc
    else:
        mols = Nv[0]
        V = Vnorm
except:
    mols = Nv[0]
    V = Vnorm

# Theta
theta = np.array(mols/mols[0])

# Capacidades caloríficas (cal/kmol)
Cpx = np.array([40, 8.38, 18])

# Entalpia padrão de reação (kcal/kmol)
dHrx = -5.9e5

# NiCpi (kcal/K)
NiCpi = np.dot(mols,Cpx)

# Temp de referência (K)
Ta = 298

# UA (kcal/min Cº)
UA = 35.85 

# Concentrações
C = np.array(mols/V)

# Definir caso de refrigeração
# PADRÃO: Refrigeração falha e não há disco
try:
    refr_case = callbacks[argv[1]]
except:
    refr_case = callbacks["rfnd"]

# Solver
sol = solve_ivp(edo,(0,120),(0,448),max_step=1)

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

file = open('resultados','w')
file.write("    t       |       X      |       T      |      Qrem    |    Qgen    |\n")
for t, X, T, Qr, Qg in zip(sol.t, sol.y[0], sol.y[1], Qrem(sol.y[1]), Qgen(sol.y[1],sol.y[0])):
    file.write("-"*60)
    file.write("\n")
    file.write(f"{t:^10.5}  |  {X:^10.5}  |  {T:^10.5}  |  {Qr:^10.5}  |  {Qg:^10.5}\n")
file.close()
