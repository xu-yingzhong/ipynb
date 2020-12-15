import matplotlib.pyplot as plt
import matplotlib as mlb
import numpy as np
import toycode1 as toy

basePath = '/srv/cosmdatb/erebos/peng/sims.TNG/TNG300-1/output'
subhaloID=0
parData = toy.SDseparateCC(basePath,snapNum=99,
                              subhaloID=subhaloID,a=1,H_z_0=1)
''' make grid '''

Sx = parData['spheroid']['Coordinates'][:,0]
Sy = parData['spheroid']['Coordinates'][:,1]
Sz = parData['spheroid']['Coordinates'][:,2]
Sm = parData['spheroid']['Masses']
Sv = parData['spheroid']['Velocities'][:,1]

Dx = parData['disc']['Coordinates'][:,0]
Dy = parData['disc']['Coordinates'][:,1]
Dz = parData['disc']['Coordinates'][:,2]
Dm = parData['disc']['Masses']
Dv = parData['disc']['Velocities'][:,1]


# get S/T
SphMass = 0
DiscMass = 0
for i in np.arange(len(parData['spheroid']['Masses'])):
	SphMass = SphMass + parData['spheroid']['Masses'][i]

for i in np.arange(len(parData['disc']['Masses'])):
	DiscMass = DiscMass + parData['disc']['Masses'][i]

ratioS_T = SphMass/(SphMass+DiscMass)

#plot
#fig1

fig1, (ax1,ax2) = plt.subplots(1,2,sharex=True,sharey=True)
fig1.suptitle('subhaloID=%d (S/T=%5.3f)' %(subhaloID,ratioS_T))

ax1.scatter(parData['spheroid']['cirPar'],parData['spheroid']['radius'],label='sph',s=10)
ax2.scatter(parData['disc']['cirPar'],parData['disc']['radius'],label='disc',s=10)
ax1.set_xlabel('$\epsilon$')
ax1.set_ylabel('$R(kpc/h)$')
ax2.set_xlabel('$\epsilon$')
ax1.legend()
ax2.legend()
'''
#fig2

fig2, (ax3,ax4) = plt.subplots(1,2,sharex=True,sharey=True)
fig2.suptitle('subhaloID=%d (edge-on)'%subhaloID)
SPolyCollec = ax3.hexbin(Sx,Sz)
DPolyCollec = ax4.hexbin(Dx,Dz)
Sl = SPolyCollec.get_offsets()[1,1]-SPolyCollec.get_offsets()[0,1]
Dl = DPolyCollec.get_offsets()[1,1]-DPolyCollec.get_offsets()[0,1]
#############################
#print(SPolyCollec.get_offsets()[0,:])
#print(SPolyCollec.get_offsets()[1,:])
#print(DPolyCollec.get_offsets()[0,:])
#print(DPolyCollec.get_offsets()[1,:])
#############################
SHexArea = np.sqrt(3)/2*Sl**2
DHexArea = np.sqrt(3)/2*Dl**2
Sm=Sm/SHexArea
Dm=Dm/DHexArea
ax3.cla()
ax4.cla()
pcm=ax3.hexbin(Sx,Sz,C=Sm,bins='log',cmap='inferno',reduce_C_function=np.sum)
ax4.hexbin(Dx,Dz,C=Dm,bins='log',cmap='inferno',reduce_C_function=np.sum)
cb = fig2.colorbar(pcm,ax=[ax3,ax4],shrink=0.6,location='bottom')
cb.set_label('$\log\Sigma\ (M_{\odot} h/kpc^2)$')
ax3.set_xlabel('$x(kpc/h)$')
ax3.set_ylabel('$y(kpc/h)$')
ax4.set_xlabel('$x(kpc/h)$')
ax3.set_title('spheroid')
ax4.set_title('disc')

#fig3

fig3, (ax5,ax6) = plt.subplots(1,2,sharex=True,sharey=True)
fig3.suptitle('subhaloID=%d (edge-on)'%subhaloID)

cmap2 = plt.cm.get_cmap('seismic')
pcm=ax5.hexbin(Sx,Sz,C=Sv,cmap=cmap2)
ax6.hexbin(Dx,Dz,C=Dv,cmap=cmap2)
cb = fig3.colorbar(pcm,ax=[ax5,ax6],shrink=0.6,spacing='proportional',location='bottom')
cb.set_label('$V\ (km/s)$')
ax5.set_xlabel('$x(kpc/h)$')
ax5.set_ylabel('$z(kpc/h)$')
ax6.set_xlabel('$x(kpc/h)$')
ax5.set_title('spheroid')
ax6.set_title('disc')
'''
plt.show()

