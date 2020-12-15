import toycode1 as toy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mlb

''' select subhalo and write in a file'''
result1 = {}
result2 = {}
Fields = ['ID','ratioS_T','sSFR']
basePath = '/srv/cosmdatb/erebos/peng/sims.TNG/TNG300-1/output'
result = toy.selectSubhalo(basePath,snapNum=99,limStarMass=10**11,
                           sSfrLog=-11,Str=1,a=1,H_z_0=1,trigger=False)
idQuenchedDisc = np.intersect1d(np.where(result['sSFR']<-11),np.where(result['ratioS_T']<0.4))

for field in Fields:
        result1[field] = result[field][idQuenchedDisc]

OderIndex = np.argsort(result1['ratioS_T'])
for field in Fields:                                #sort the result1 according to S/T
        result1[field] = result1[field][OderIndex]

for field in Fields:
        result2[field] = np.delete(result[field],idQuenchedDisc)
'''archive'''
for i in np.arange(len(result1['ID'])):
    f = open('/srv/cosmdatb/erebos/yingzhong/selData/data/result_04','a') #sometimes change 'a' to 'w'
    f.write('\nID: %d ratioS_T: %5.3f sSFR: %f' % (result1['ID'][i],result1['ratioS_T'][i],result1['sSFR'][i]))
    f.close()

np.savez('/srv/cosmdatb/erebos/yingzhong/selData/data/result_04.npz',
         ID = result1['ID'],ratioS_T=result1['ratioS_T'],sSFR=result1['sSFR'])

np.savez('/srv/cosmdatb/erebos/yingzhong/selData/data/result_full.npz',
         ID = result['ID'],ratioS_T=result['ratioS_T'],sSFR=result['sSFR'])

'''plot SFR-S/T '''

fig1,ax1= plt.subplots(1,1)
#fig1.suptitle('S/T < 0.5')

ax1.scatter(x=result2['ratioS_T'],y=result2['sSFR'],s = 10)
ax1.scatter(x=result1['ratioS_T'],y=result1['sSFR'],s = 10,label='$sSFR < 10^{-11}\ yr^{-1}$ & $S/T<0.4$')
#ax1.axhline(y=10**(-11),c='r')
#ax2.hist(x=result1['ratioS_T'],bins='auto')
#cb = fig1.colorbar(pcm[3],ax=ax1,shrink=0.6,spacing='proportional')
#ax1.set_ylim([-0.05,1])
ax1.set_xlabel('$S/T$')
ax1.set_ylabel('$\log_{10}sSFR $')
#ax2.set_xlabel('$S/T$')
#ax2.set_ylabel('$N$')

ax1.set_title('NO limitation on sSFR')
#ax2.set_title('$sSFR < 10^{-11}\ yr^{-1}$')

ax1.legend()
plt.savefig('/srv/cosmdatb/erebos/yingzhong/selData/fig/fig')


      
