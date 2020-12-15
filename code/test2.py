import numpy as np
import matplotlib.pyplot as plt

fMmass = np.load('/srv/cosmdatb/erebos/yingzhong/selData/data/haveBH.npz')
f1     = np.load('/srv/cosmdatb/erebos/yingzhong/selData/data/haveBH_025.npz')

a1 = 8.46
b1 = 1.05

a2 = 8.69
b2 = 1.16

x   = np.linspace(9.5,12,20)
y1  = b1*x + a1 - 11*b1
y2  = b2*x + a2 - 11*b2

fig, ax = plt.subplots()
ax.scatter(np.log10(fMmass['M_bulge']),np.log10(fMmass['BH_CenMass']),s=1)
ax.scatter(np.log10(f1['M_bulge']),np.log10(f1['BH_CenMass']),s=10,
           label='$sSFR < 10^{-11}\ yr^{-1}$ & $S/T<0.25$')
ax.plot(x,y1,'-k',label = 'McConnell&Ma 2013')
ax.plot(x,y2,'--k',label = 'Kormendy&Ho 2013')

ax.set_xlabel('$\log_{10}(M_{bulge}/M_{\odot})$')
ax.set_ylabel('$\log_{10}(M_{BH}/M_{\odot})$')
ax.legend()

plt.show()
