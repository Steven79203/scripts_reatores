#!/usr/bin/python

# Dados do problema 11.3a - CSTR
# Reação elementar: A + B -> C

# a) Compare com dois reatores CSTR em série

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# == FUNÇÕES AUXILIARES == # 

def rx(X,T):
    return -(ki(k,T) * Ca0 * (1-X) * Ca0 * (theta_b - X))

def ki(k,T):
    return k * np.exp(E/R * (1/Tr - 1/T))

def dHrx(T):
    return deltaHrx + dCp * (T - Tr)

def theta_cpi(X):
    Ca = Ca0*(1-X)
    Cb = Ca0*(theta_b-X)
    Cc = Ca0*X
    return Cpa + Cpb*Cb/Ca + Cpc*Cc/Ca

# == FUNÇÃO OBJETIVO == #

def func(x):
    X = x[0]
    T = x[1]
    F = V*(-rx(X,T)) - Fa0*X
    G = -theta_cpi(X)*(T-T0) - dHrx(T)*X

    return [F,G]

# == DADOS == # 

a = 1
b = 1
c = 1

V = 250 # dm³

Tr = 273

v0  = 2   # dm³/s
Ca0 = 0.1 # mol/dm³
Fa0 = v0 * Ca0
Fb0 = Fa0

theta_b = Fb0/Fa0

Ha = -20000 # cal/mol
Hb = -15000 # cal/mol
Hc = -41000 # cal/mol
deltaHrx = (c/a)*Hc - (b/a)*Hb - Ha

Cpa = 15 # cal/mol.K
Cpb = Cpa 
Cpc = 30 # cal/mol.K
dCp = (c/a) * Cpc - (b/a)*Cpb - Cpa

k = 0.01 # dm³/mol.s
E = 10000 # cal/mol
R = 1.987 # cal/mol.K

#UA  = 20 # cal/m³.s.K
#mc  = 50 # g/s 
#Ta0 = 450 # K
#Cpc = 1 # cal/g.K
 
T0 = 273+27
x0 = [0.9, 350]
res = fsolve(func, x0)
print(res)

# Recalcular dados
Ca0 = Ca0*(1-res[0])
Fa0 = Ca0*v0
Fb0 = Fa0
theta_b = Fb0/Fa0
T0 = res[1]

x0 = res
res = fsolve(func,x0)
print(res)

