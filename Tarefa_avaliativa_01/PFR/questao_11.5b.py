#!/usr/bin/python 

# Problema 11.5b

# Reação endotérmica A -> B +C 

# V = 500 dm³ 

# == FUNÇÃO AUXILIAR == #
def ki(T):
    return np.exp(34.34 - 34.222/T)

def rx(X,T):
    Ca = Ca01*(1-X)
    return -(ki(T)*Ca)

def dHrx(T):
    return deltaHrx + dCp*(T-Tr)

def pfr(t,F):
    X = F[0]
    T = F[1]
    dXdV = -rx(X,T)/Fa0
    dTdV = (-rx(X,T)*(-dHrx(T)))/(Fa0 * (theta_cpi + X*dCp))
    return [dXdV, dTdV]

# == DADOS DO PROBLEMA == # 

R = 0.082 # atm.L/mol.K
T = 1100  # K

Fa0  = 10 # mol/min
Pa0  = 2  # atm
Ca0 = P/(R*T)

Ci0  = P/(R*T)

Ct   = P/(R*T)
Ca01 = 0.5*Ct

theta_i = Ci0/Ca0

Cpa = 170 #J/mol.k
Cpb = 90  #J/mol.K
Cpc = 80  #J/mol.K
dCp = Cpc + Cpb - Cpa

deltaHrx = 80000 #J/mol

Tr = 273  #K 
T0 = 1100 #K
Ta = 0    #L

res = solve_ivp(pfr,(0,500),(0,T0))
