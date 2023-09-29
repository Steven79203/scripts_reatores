#!/usr/bin/env python 

# DEQ1032 - Engenharia das reações químicas avançadas
# Reator não-adiabático com trocador de calor 
# Aula 02 - 31 de Março de 2023


from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt 
import numpy as np

# == Funções auxiliares == #
# Parâmetro cinético k 
def kin(T):
    return 4.71 * 10**9 * np.exp(-18000/ (R*T))

# Concentração em X
def C(C0, X):
    return C0 * (1-X)

# Taxa de reação
def r(k,C):
    return -k * C

# delta Hrx (Entalpia de reação)
def dHrx(T):
    return -20202 + dCp * (T - TR)

# Calor removido 
def Qrem(T):
    return mc * Cpfluido * (T - Ta) * (1 - np.exp(-UA/mc/Cpfluido))

# Calor gerado
def Qgen(T,X):
    return -dHrx(T) * ( -r(kin(T),C(Ca0,X)) * V )

# == Definição da EDO == #
def edo(t, F):

    Qremi = Qrem(F[1])
    Qgeni = Qgen(F[1], F[0])

    dXdt = -r(kin(F[1]),C(Ca0,F[0])) / Ca0
    dTdt = (Qgeni - Qremi)/(Na0 * (ThetaCp + dCp * F[0]))

    # Retorno do sistema de EDOs
    return [dXdt,dTdt]

# Mols dos componentes
Na0 = 54.8
Nb0 = 555 
Nc0 = 0
Ni0 = 98.8

# 1 - Capacidade calorífica (cal/mol*K)
Cpa = 35
Cpb = 18
Cpc = 46
Cpi = 19.5

# Volume (dm³)
V = 18

# 2 - Entalpia de formação (cal/mol)
# Desconsiderar caso dHrx tenha sido fornecido
Ha = -12667.3
Hb = -68300
Hc = -119502.86  
#Ha = -66000
#Hb = -123000
#Hc = -226000

TR = 293 #K
Ta = 285 #K

# Reatores Batelada adiabático 
T0 = 286
theta_a = Na0/Na0
theta_b = Nb0/Na0
theta_c = Nc0/Na0
theta_i = Ni0/Na0

# Balanço molar 
R = 1.987 #cal/mol.K
a = 1
b = 1
c = 1
Ca0 = Na0/V

# Balanço de energia 
dH0 = ((c/a)*Hc) - ((b/a)*Hb) - (Ha) 
dCp = (c/a) * Cpc - (b/a) * Cpb - Cpa
ThetaCp = (theta_a * Cpa) + (theta_b * Cpb) + (theta_c * Cpc) + (theta_i * Cpi)

# Coeficiente global de transferência de calor 
UA = 10 # cal/
mc = 10         #g/s
Cpfluido = 4.16 #cal/g*K

# Chamada da função solve_ivp
sol = solve_ivp(edo,(0,4000),(0,T0),max_step=1)

# Plot do gráfico
#Conversão
fig1, ax1 = plt.subplots()
ax1.plot(sol.t,sol.y[0],'g-',label='Conversão')
ax1.set_xlabel('tempo (seg)')
ax1.set_ylabel('X (mol/mol)')
ax1.set_xlim([0,4000])
ax1.grid()
ax1.legend()

#Temperatura
fig2, ax2 = plt.subplots()
ax2.plot(sol.t,sol.y[1],'b-',label='Temperatura')
ax2.set_xlabel('tempo(seg)')
ax2.set_ylabel('T (K)')
ax2.set_xlim([0,4000])
ax2.set_ylim([273,400])
ax2.grid()
ax2.legend()

#Calor gerado e removido
fig3,ax3 = plt.subplots()
ax3.plot(sol.t, Qgen(sol.y[1],sol.y[0]),color="m",label="Calor gerado")
ax3.plot(sol.t, Qrem(sol.y[1]),color="b",label='Calor removido')
ax3.set_ylabel("Calor (cal)")
ax3.set_xlabel("tempo (seg)")
ax3.set_xlim([0,4000])
ax3.grid()
ax3.legend()

plt.show()

file = open('resultados.txt','w')
file.write("    t     |       X      |      T     \n")
for t, X, T in zip(sol.t, sol.y[0], sol.y[1]):
    file.write("-"*40)
    file.write("\n")
    file.write(f"{t:^10.5}  |  {X:^10.5}  |  {T:^10.5}\n")
file.close()
