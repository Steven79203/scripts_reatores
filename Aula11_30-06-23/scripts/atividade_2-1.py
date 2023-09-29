#!/usr/bin/python

from scipy.integrate import solve_ivp
import matplotlib.pyplot as plt

# EDO
def edo(t, F):
    Ca, Cb, Cc, Cd, Ce, Ch, Cf, Ck = F[::]
    
    return [dCadt(Ca,Cb),
            dCbdt(Ca,Cb),
            dCcdt(Cd),
            dCddt(Ca,Cb,Cd),
            dCedt(Cd,Ck,Ch,Ce,Cf),
            dChdt(Ck,Ch,Ce,Cf),
            dCfdt(Ck,Ch,Ce,Cf),
            dCkdt(Ce,Cf,Ck,Ch)]

# EDOs
dCadt = lambda Ca,Cb         : -r1(Ca,Cb)
dCbdt = lambda Ca,Cb         : -r1(Ca,Cb)
dCcdt = lambda Cd            : 2*r2(Cd)
dCddt = lambda Ca,Cb,Cd      : (r1(Ca,Cb)-r2(Cd))
dCedt = lambda Cd,Ck,Ch,Ce,Cf: r2(Cd) + r6(Ck,Ch) - r3(Ce,Cf)
dChdt = lambda Ck,Ch,Ce,Cf   : -2*r6(Ck,Ch) + 2*r3(Ce,Cf)
dCfdt = lambda Ck,Ch,Ce,Cf   : r5(Ck) + r4(Ck) + r6(Ck,Ch) - r3(Ce,Cf)
dCkdt = lambda Ce,Cf,Ck,Ch   : r3(Ce,Cf) - r4(Ck) -r5(Ck) - r6(Ck,Ch)

# Auxiliares
r1 = lambda Ca,Cb: k1*Ca*Cb
r2 = lambda Cd   : k2*Cd
r3 = lambda Ce,Cf: k3*Ce*Cf
r4 = lambda Ck   : k4*Ck
r5 = lambda Ck   : k5*Ck
r6 = lambda Ck,Ch: k6*Ck*Ch**2
rEMIT = lambda Ck: r4(Ck)*QY

# Constantes
k1 = 1.485e7  #dm³/mol/s
k2 = 1.485e7  #s-1
k3 = 891      #dm³/mol/s
k4 = 0.05     # s^-1
k5 = 0.111    # s^-1
k6 = 1.782e12 # dm³/mol/s
QY = 200000

Ci0 = [0.7, 0.8, 0, 0, 0, 0, 0.0005, 0]
tf = 3600*25*5

res = solve_ivp(edo, (0,tf), Ci0, method="Radau")

print(res)

plt.plot(res.t, rEMIT(res.y[-1]))
plt.title('Luminosidade x tempo')
plt.xlabel('tempo(s)')
plt.ylabel('Emissão')
plt.grid()
plt.show()
