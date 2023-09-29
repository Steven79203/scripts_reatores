#!/usr/bin/env python 

# DEQ1032 - Engenharia das reações químicas avançadas
# Aula 01

from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt 
import numpy as np

# Definição da EDO exigida para o solve_ivp
def edo(t, F):
    """
    Todas as equações que dependem das variáveis T e X 
    precisam estar dentro da definição da EDO
    O "Vetor" F referem-se às funções T e X, sendo:
        F[0] = X
        F[1] = T
    """
    
    # Parâmetros cinéticos
    #k = 16.71 * 10**12 * np.exp(-18000 / (R*F[1]))
    k = 4.71 * 10**9 * np.exp(-18000 / (R*F[1]))

    # Concentração
    Ca = Ca0 * (1-F[0])

    # Taxa de reação
    ra = -k * Ca

    # deltaHrx
    #deltaHrx = dH0 + dCp * (F[1]-TR)
    deltaHrx = -20202 + dCp * (F[1]-TR)

    # EDOs
    dXdt = -ra / Ca0
    dTdt = ((-deltaHrx) * (-ra))/(Ca0 * (ThetaCp + dCp * F[0]))

    # Retorno do sistema de EDOs
    return [dXdt,dTdt]

# Mols dos componentes
Na = 54.8
Nb = 555 
Nc = 0
Ni = 98.8

# 1 - Capacidade calorífica (cal/mol*K)
Cpa = 35
Cpb = 18
Cpc = 46
Cpi = 19.5

# Volume (dm³)
V = 4

# 2 - Entalpia de formação
Ha = -12667.3
Hb = -68300
Hc = -119502.86
#Ha = -66000
#Hb = -123000
#Hc = -226000

TR = 293 #K

# Reatores Batelada adiabático 
T0 = 286
theta_a = Na/Na
theta_b = Nb/Na
theta_c = Nc/Na
theta_i = Ni/Na

# Balanço molar 
R = 1.987 #cal/mol.K
a = 1
b = 1
c = 1
Ca0 = 0.13

# Balanço de energia 
dH0 = ((c/a)*Hc) - ((b/a)*Hb) - (Ha) 
dCp = (c/a) * Cpc - (b/a) * Cpb - Cpa
ThetaCp = (theta_a * Cpa) + (theta_b * Cpb) + (theta_c * Cpc) + (theta_i * Cpi)

# Chamada da função solve_ivp
sol = solve_ivp(edo,(0,4000),(0,T0))


# Plot do gráfico
fig1, ax1 = plt.subplots()
ax1.plot(sol.t,sol.y[0],'g-',label='Conversão')
ax1.legend()

fig2, ax2 = plt.subplots()
ax2.plot(sol.t,sol.y[1],'b-',label='Temperatura')
ax2.legend()

plt.show()

# 
print('t   -   X   -   T')
for t, X, T in zip(sol.t[::50], sol.y[0][::50], sol.y[1][::50]):
    print(f'{t:12.10}  -  {X:12.10}  -  {T:12.10}')
