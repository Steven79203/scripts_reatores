#!/usr/bin/env python

from scipy.integrate import trapezoid
import numpy as np
import matplotlib.pyplot as plt

t = np.array([
0,
31.92,
39.17,
45.17,
51.44,
57.59,
64.29,
70.44,
76.63,
82.61,
89.22,
95.25,
101.9,
109.19,
116.54,
123.72,
138.31])

C = np.array([
0,
0.8592298,
7.3365006,
4.4343468,
2.7218958,
1.952795 ,
1.2858404,
0.9493588,
0.4927052,
0.5287568,
0.4025762,
0.330473 ,
0.2523612,
0.1922752,
0.1381978,
0.1021462,
0.0360516,
])

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

print("===========")
print("E(t): ",res)
print("tm", tm)
print("sig²", sig_2)

print("")
print(f"{E}\n{Et}\n{Etm2}")

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



