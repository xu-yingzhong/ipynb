import matplotlib.pyplot as plt
import matplotlib as mlb
import numpy as np
import toycode1 as toy
import seaborn as sns


'''============================================================================'''
def massDis(ax,x,y,m):
    PolyCollec = ax.hexbin(x,y)

    l = PolyCollec.get_offsets()[1,1]-PolyCollec.get_offsets()[0,1]

    HexArea = np.sqrt(3)/2*l**2
    m=m/HexArea
    ax.cla()

    pcm=ax.hexbin(x,y,C=m,cmap='inferno',reduce_C_function=np.sum)

    return pcm
    
'''=========================================================================='''

'''==========================================================================='''
def velDis(ax,x,y,v_z):
    cmap2 = plt.cm.get_cmap('seismic')
    pcm=ax.hexbin(x,y,C=v_z,cmap=cmap2)

    return pcm

'''============================================================================'''
sns.set(color_codes=True)

basePath = '/srv/cosmdatb/erebos/peng/sims.TNG/TNG300-1/output'
subhaloID = 0
snapNum = 99
#Type = ['spheroid','disc']
#string = Type[1]
#parData = toy.SDseparateCC(basePath,snapNum=135,subhaloID=subhaloID,a=1,H_z_0=1)['%s' %string]

parData = toy.cirparameter(basePath,snapNum=snapNum,ID=subhaloID,a=1,H_z_0=1)
''' make grid '''

c = parData['Coordinates']
m = parData['Masses']
v = parData['Velocities']
#cirPar = parData['cirPar']

''' get B/T'''

Nr = 2
Nc = 3
indexX = [0,0,1]
indexY = [1,2,2]
indexZ = [2,1,0]
stringX = ['x','x','y']
stringY = ['y','z','z']

fig,axs = plt.subplots(Nr,Nc,sharex=True,sharey=True)
#fig.suptitle('0.1 ID: %d (%s)' %(subhaloID,string))
fig.suptitle('ID: %d' %(subhaloID))
for i in np.arange(Nr):
    pcm = [0,0,0]
    for j in np.arange(Nc):
        ax = axs[i,j]
        if i == 0:
            x = c[:,indexX[j]]
            y = c[:,indexY[j]]
            pcm[j] = massDis(ax,x,y,m)
        else:
            x   = c[:,indexX[j]]
            y   = c[:,indexY[j]]
            v_z = v[:,indexZ[j]]
            pcm[j] = velDis(ax,x,y,v_z)
            
        ax.set_xlabel(stringX[j]+'$\ (kpc/h)$')
        ax.set_ylabel(stringY[j]+'$\ (kpc/h)$')

    vmin = min(image.get_array().min() for image in pcm)
    vmax = max(image.get_array().max() for image in pcm)
    if i == 0:
        norm = mlb.colors.LogNorm(vmin=vmin, vmax=vmax)
    else:
        if abs(vmin) < vmax:
            vmin = -1*vmax
        else:
            vmax = abs(vmin)
        norm = mlb.colors.Normalize(vmin=vmin, vmax=vmax)
    for im in pcm:
        im.set_norm(norm)
    cb = fig.colorbar(pcm[0],ax=[axs[i,:]])
    if i == 0:
        cb.set_label('$\log\Sigma\ (M_{\odot} h/kpc^2)$')
    else:
        cb.set_label('$V\ (km/s)$')
                
plt.savefig('/z/yingzhong/fig/fig')
#plt.show()
