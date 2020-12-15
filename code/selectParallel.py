import illustris_python1 as il
import toycode1 as toy
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mlb
from multiprocessing import Pool

basePath = '/srv/cosmdatb/erebos/peng/sims.TNG/TNG300-1/output'
snapNum = 99
h = 0.677
limStarMass = 10**10.5

def selectFunc(ID):
        a = 1
        H_z_0 = 1
        
        return toy.selectSubhalo(basePath,snapNum,a,H_z_0,ID,limStarMass)

#def changeFormat(i):
#    idarray  = np.append(idarray,results[i][0],axis=0)
#    sSFRarray = np.append(sSFRarray,results[i][2],axis=0)
#    S_T_arr  = np.append(S_T_arr,results[i][1],axis=0)
    
''' select subhalo and write in a file'''
result = {}
result1 = {}
result2 = {}
idarray = np.zeros(0,dtype=int)
sSFRarray = np.zeros(0)
S_T_arr = np.zeros(0)
Fields = ['ID','ratioS_T','sSFR']

roughLim = (limStarMass*h) / 10**10
firstSubId = il.groupcat.loadHalos(basePath,snapNum,fields='GroupFirstSub')
firstSubId1 = firstSubId[np.where(firstSubId != -1)]
SubhaloSMass = il.groupcat.loadSubhalos(basePath,snapNum,fields='SubhaloMassType')[:,4]
firSubhaloSMass = SubhaloSMass[firstSubId1]
firstSubId = firstSubId1[np.where(firSubhaloSMass >= roughLim)]

p = Pool(40)
results = p.map(selectFunc,firstSubId)
p.close()
p.join()
        
for i in np.arange(len(results)):
    idarray  = np.append(idarray,results[i][0],axis=0)
    sSFRarray = np.append(sSFRarray,results[i][2],axis=0)
    S_T_arr  = np.append(S_T_arr,results[i][1],axis=0)

result['ratioS_T'] = S_T_arr
result['ID'] = idarray
result['sSFR'] = sSFRarray
    
idQuenchedDisc = np.intersect1d(np.where(result['sSFR']<-11),np.where(result['ratioS_T']<0.4))

for field in Fields:
        result1[field] = result[field][idQuenchedDisc]

OderIndex = np.argsort(result1['ratioS_T'])
for field in Fields:                                #sort the result1 according to S/T
        result1[field] = result1[field][OderIndex]

for field in Fields:
        result2[field] = np.delete(result[field],idQuenchedDisc)

'''archive'''

with open('/srv/cosmdatb/erebos/yingzhong/selData/data/result_04P','a') as f:
    for i in np.arange(len(result1['ID'])):
        f.write('\nID: %d ratioS_T: %5.3f sSFR: %f' % (result1['ID'][i],
                                                       result1['ratioS_T'][i],result1['sSFR'][i]))

np.savez('/srv/cosmdatb/erebos/yingzhong/selData/data/result_04P.npz',
         ID = result1['ID'],ratioS_T=result1['ratioS_T'],sSFR=result1['sSFR'])

np.savez('/srv/cosmdatb/erebos/yingzhong/selData/data/result_full_P.npz',
         ID = result['ID'],ratioS_T=result['ratioS_T'],sSFR=result['sSFR'])




'''plot SFR-S/T '''
'''
fig1,ax1= plt.subplots(1,1)
#fig1.suptitle('S/T < 0.5')

ax1.scatter(x=result2['ratioS_T'],y=result2['sSFR'],s = 10)
ax1.scatter(x=result1['ratioS_T'],y=result1['sSFR'],s = 10,label='$sSFR < 10^{-11}\ yr^{-1}$ & $S/T<0.25$')
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
'''

      
