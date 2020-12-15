import numpy as np
import matplotlib.pyplot as plt

###########################################################
'''plot M_bulge versus M_bh'''
###########################################################

filedsMass = ['BH_CenMass', 'M_star', 'BH_ID', 'subHaloID', 'M_dm_200c']
filedsMassFull = ['BH_CenMass', 'M_star', 'BH_ID', 'subHaloID', 'M_dm_200c','M_bulge']
fieldsRatio = ['sSFR', 'ratioS_T', 'ID']
f = {}
f1 = {}

'''load data'''
fMmass = np.load('/srv/cosmdatb/erebos/yingzhong/selData/data/resultMass.npz')
fRatio = np.load('/srv/cosmdatb/erebos/yingzhong/selData/data/result_full_P.npz')
f_025_Ratio = np.load('/srv/cosmdatb/erebos/yingzhong/selData/data/result_04P.npz')

fMmass = dict(fMmass)
fRatio = dict(fRatio)
f_025_Ratio = dict(f_025_Ratio)


'''delete subhalo without bh'''
nullBH_IndexMass = np.where(fMmass['BH_ID']==-1)
nullSubHaloId = fMmass['subHaloID'][nullBH_IndexMass]

for field in filedsMass:
    fMmass[field] = np.delete(fMmass[field],nullBH_IndexMass)

for subHaloId in nullSubHaloId:
    for field in fieldsRatio:
        fRatio[field] = np.delete(fRatio[field],np.where(fRatio['ID']==subHaloId))
    
'''sort'''
for field in fieldsRatio:
    f[field] = np.zeros(0)

for SubhaloId in fMmass['subHaloID']:
    for field in fieldsRatio:
        if field == 'ratioS_T' and fRatio[field][np.where(fRatio['ID']==SubhaloId)]>1:
            f[field] = np.append(f[field],1)
        else:
            f[field] = np.append(f[field],fRatio[field][np.where(fRatio['ID']==SubhaloId)])
            
f['ID'] = f['ID'].astype(int)
fMmass['subHaloID'] = fMmass['subHaloID'].astype(int)

'''calculate M_Bulge'''
index = np.zeros(0,dtype=int)
fMmass['M_bulge'] = f['ratioS_T']*fMmass['M_star']

for SubhaloId in f_025_Ratio['ID']:
    index = np.append(index,np.where(fMmass['subHaloID']==SubhaloId)[0])

for field in filedsMassFull:
    f1[field] = fMmass[field][index]

for field in fieldsRatio:
    if field != 'ID':
        f1[field] = f[field][index]
        

'''archive'''
np.savez('/srv/cosmdatb/erebos/yingzhong/selData/data/haveBH.npz',
         BH_CenMass=fMmass['BH_CenMass'],M_star=fMmass['M_star'],
         BH_ID=fMmass['BH_ID'],subHaloID=fMmass['subHaloID'],
         M_dm_200c=fMmass['M_dm_200c'],M_bulge=fMmass['M_bulge'],
         sSFR=f['sSFR'],ratioS_T=f['ratioS_T'])                         # all S/T <1

np.savez('/srv/cosmdatb/erebos/yingzhong/selData/data/haveBH_04.npz',
         BH_CenMass=f1['BH_CenMass'],M_star=f1['M_star'],
         BH_ID=f1['BH_ID'],subHaloID=f1['subHaloID'],
         M_dm_200c=f1['M_dm_200c'],M_bulge=f1['M_bulge'],
         sSFR=f1['sSFR'],ratioS_T=f1['ratioS_T']) 

'''plot'''
'''
fig, ax = plt.subplots()
ax.scatter(np.log10(fMmass['M_bulge']),np.log10(fMmass['BH_CenMass']))
ax.scatter(np.log10(f1['M_bulge']),np.log10(f1['BH_CenMass']))
ax.set_xlabel('M_bulge')
ax.set_ylabel('BH_CenMass')

plt.show()
'''