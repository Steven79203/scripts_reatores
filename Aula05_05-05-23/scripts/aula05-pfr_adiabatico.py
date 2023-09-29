#!/usr/bin/python 

# Reação 
# A + B -> C

from scipy.integrate import solve_ivp 
import matplotlib.pyplot as plt
import numpy as np

def Temp(X):
    N = -X*deltaHrx + s_theta_cpi*T0 + X*deltaCp*Tr
    D = s_theta_cpi + X*deltaCp
    return N/D

def Xe(T):
    KC = Kc(Kc2, T)
    return KC / (1+KC)

def Kc(K,T):
    return K * np.exp(deltaHrx/R * (1/T2 - 1/T))

def ki(k,T):
    return k * np.exp(E/R * (1/T1 - 1/T))

def rt(T,X):
    Ca = Ca0*(1-X)
    Cb = Ca0*(theta_b + (b/a)*X)
    return -ki(k1,T) * (Ca - Cb/Kc(Kc2,T))

def edo(t,F):
    X = F[0]
    dXdV = -rt(Temp(X),X)/Fa0
    return dXdV


# = Dados = # 
Ft0 = 163.3    # kmol_a/h
ya0 = 0.9      # kmol_a/kmol
yi0 = 0.1      # kmol_a/kmol 
Fa0 = ya0*Ft0  # kmol_a/h
Ca0 = 9.3      # kmol_a/dm³

# Coeficientes estequimétricos
a = 1
b = 1
c = 1

# Capacidades caloríficas (J/mol.K)
Cpa = 141
Cpb = 141
Cpi = 161
deltaCp = (b/a)*Cpb - Cpa

# Theta
theta_b = 0
theta_i = 0.1/0.9

# Somatheta
s_theta_cpi = Cpa + theta_i*Cpi + theta_b*Cpb

# Temperaturas (K)
T0 = 330 
T1 = 360
T2 = 333
Tr = 330

# Dados cinéticos e de equilíbrio
k1  = 31.1    # h^-1
E   = 65.7e3  # J/mol
R   = 8.314   # J/mol.K
Kc2 = 3.03

deltaHrx = -6900 # J/mol

# Volumes
V0 = 0
Vf = 5

res = solve_ivp(edo,(0,Vf),(0,0))

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()

t = res.t
X = res.y[0]

ax1.plot(t, X)              # Conversão
ax1.plot(t, Xe(Temp(X)))    # Conversão de equilíbrio
ax1.set_title('Conversão')

ax2.plot(t, -rt(Temp(X),X))  # Taxa de reação
ax2.set_title('Taxa de reação (-ra)')

ax3.plot(t, Temp(X))        # Temperatura
ax3.set_title('Temperatura(K)')

plt.show()
