#!/usr/bin/env python

from scipy.integrate import trapezoid
import numpy as np
import matplotlib.pyplot as plt

t = np.array([0, 18.49, 21.65, 25.18, 28.48, 32.01, 34.37, 37.14, 40.81, 43.83, 46.83, 50.34, 53.78, 57.79, 61.64])

C = np.array([0, 2.11, 7.63, 5.66, 2.17, 1.03, 0.77, 0.65, 0.59, 0.54, 0.51, 0.47, 0.43, 0.40, 0.36])

k = 0.1

# Área da curva C x t
res = trapezoid(C,t)

# E(t)
E = C/res

# E(t) * t 
Et = E*t

# Tempo médio
tm = trapezoid(Et,t)

# E(t-tm)² 
Etm2 = E*(t-tm)**2

# sigma**2
sig_2 = trapezoid(Etm2,t)

print("E(t): ",res)
print("tm", tm)
print("sig²", sig_2)
exit()

# Conversão (X)
X = 1 - np.exp(-k*t)

# E(t)*X(t)
EX = E*X

# X_medio
Xm = trapezoid(EX,t)

fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()
fig4, ax4 = plt.subplots()


ax1.plot(t,C)
ax1.set_title("C x t")

ax2.plot(t,E)
ax2.set_title("E(t) x t")

ax3.plot(t,Et)
ax3.set_title("E(t).t x t")

ax4.plot(t,Etm2)
ax4.set_title("E(t)(t-tm)**2 x t")

plt.show()



