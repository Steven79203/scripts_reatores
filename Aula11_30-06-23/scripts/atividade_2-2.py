#!/usr/bin/python 

from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

def edo(t,F):
    Ca,Cb,Cc,Cd,Ce,Ch,Cf,Ck = F[::]
    return [dCadt(Ca,Cb),
            dCbdt(Ca,Cb),
            dCcdt(Ca,Cb,Cd),
            dCddt(Ca,Cb,Cd),
            dCedt(Cd,Cf,Ce),
            dChdt(Cf,Ce),
            dCfdt(Cf,Ce,Ck),
            dCkdt(Cf,Ce,Ck)]

# EDOs
dCadt = lambda Ca,Cb: -r1(Ca,Cb)
dCbdt = lambda Ca,Cb: -r1(Ca,Cb)
dCcdt = lambda Ca,Cb,Cd: (r1(Ca,Cb)+ r2(Cd))
dCddt = lambda Ca,Cb,Cd: (r1(Ca,Cb) - r2(Cd))
dCedt = lambda Cd,Cf,Ce: (r2(Cd) - r3(Cf,Ce))
dChdt = lambda Cf,Ce: 2*r3(Cf,Ce)
dCfdt = lambda Cf,Ce,Ck: r5(Ck) - r3(Cf,Ce)
dCkdt = lambda Cf,Ce,Ck: r3(Cf,Ce) - r5(Ck)

# Auxiliares
r1 = lambda Ca,Cb: k1*Ca*Cb
r2 = lambda Cd   : k2*Cd
r3 = lambda Cf,Ce: k3*Cf*Ce
r5 = lambda Ck   : k5*Ck

# Constantes
k1 = 1.485
k2 = 0.1485
k3 = 0.00891
k5 = 0.00111

QY = 12964
rEMIT = lambda Cf,Ce: r3(Cf,Ce)*QY

tf = 3600*24*5
Ci0 = [0.7, 0.8, 0, 0, 0, 0, 0.0005, 0]

res = solve_ivp(edo, (0,tf),Ci0,method="Radau")

plt.plot(res.t, rEMIT(res.y[6],res.y[4]))
plt.title('Luminosidade x tempo')
plt.xlabel('tempo(s)')
plt.ylabel('Emissão')
plt.grid()
plt.show()
