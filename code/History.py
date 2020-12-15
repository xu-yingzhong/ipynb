import numpy as np
import illustris_python1 as il
from multiprocessing import Pool
from functools import reduce
from scipy import integrate

def TreeAve(ID):
    
    h =  0.6774
    def AgeRedshift(x):
    
        Omega_L0 = 0.6911
        Omega_r0 = 0       #4.2 * 10**(-5) / h**2
        Omega_m0 = 0.3089

        def E(z):
            return np.sqrt(Omega_L0 - Omega_r0*(1+z)**2 + Omega_m0 * (1+z)**3 + Omega_r0*(1+z)**4)
            #return np.sqrt(Omega_L0 + Omega_m0 * (1+z)**3 + Omega_r0*(1+z)**4)
        def f(z):
            return 9.785 / (h*E(z)*(1+z))
    

        return integrate.quad(lambda z:f(z),x,np.inf)[0]

    
    SnapRedshift = dict(np.load('/srv/cosmdatb/erebos/yingzhong/selData/data/SnapRedshift.npz'))['Redshift']
    
    result     = {}
    basePath   = '/srv/cosmdatb/erebos/peng/sims.TNG/TNG300-1/output'
    snapNum    = 99
    Fields     = ['SnapNum','SubhaloMassType','SubhaloBHMass']
    Age        = np.zeros(0)
    tree       = il.sublink.loadTree(basePath, snapNum, id=ID, fields=Fields, onlyMPB=True)

    tree['SubhaloMassType'] = tree['SubhaloMassType'] * 10**10 /h  # unit: M_{\odot}
    tree['SubhaloBHMass']   = tree['SubhaloBHMass'] * 10**10 /h
    
    for i in np.arange(len(tree['SnapNum'])):
        Age = np.append(Age,AgeRedshift(SnapRedshift[tree['SnapNum']][i]))
        
    result['M_gas_max']   = np.max(tree['SubhaloMassType'][:,0])    
    index                 = np.where(tree['SubhaloMassType'][:,0] == result['M_gas_max'])[0][0] 
    result['Age_gas_max'] = tree['SnapNum'][index]
    result['M_star_max']  = tree['SubhaloMassType'][:,4][index]
    result['M_bh_max']    = tree['SubhaloBHMass'][index]
    
    result['M_gas_after']  = np.median(tree['SubhaloMassType'][0:index+1,0])
    result['M_star_after'] = np.median(tree['SubhaloMassType'][0:index+1,4])
    result['M_bh_after']   = np.median(tree['SubhaloBHMass'][0:index+1])
    
    result['M_gas_befor']  = np.median(tree['SubhaloMassType'][index:,0])
    result['M_star_befor'] = np.median(tree['SubhaloMassType'][index:,4])
    result['M_bh_befor']   = np.median(tree['SubhaloBHMass'][index:])
    
    result['PresentSubHaloID'] = int(ID)
    result['BirthTime']        = tree['SnapNum'][-1]
    
    return result


fMmass = dict(np.load('/srv/cosmdatb/erebos/yingzhong/selData/data/haveBH.npz'))
f1     = dict(np.load('/srv/cosmdatb/erebos/yingzhong/selData/data/haveBH_025.npz'))
Fields = ['M_gas_max','Age_gas_max','M_star_max','M_bh_max','M_gas_after','M_star_after',
         'M_bh_after','M_gas_befor','M_star_befor','M_bh_befor','PresentSubHaloID','BirthTime']

index  = reduce(np.intersect1d,
                (np.where(np.array(fMmass['sSFR']) > -11),
                 np.where(np.array(fMmass['ratioS_T']) <= 0.25)))

IDform   = fMmass['subHaloID'][index]
IDquench = f1['subHaloID']

res = {'M_gas_max':np.zeros(0),'Age_gas_max':np.zeros(0,dtype=int),'M_star_max':np.zeros(0),
       'M_bh_max':np.zeros(0),'M_gas_befor':np.zeros(0),'M_star_befor':np.zeros(0),
       'M_bh_befor':np.zeros(0),'M_gas_after':np.zeros(0),'M_star_after':np.zeros(0),
       'M_bh_after':np.zeros(0),'PresentSubHaloID':np.zeros(0,dtype=int),
       'BirthTime':np.zeros(0,dtype=int)}

if __name__ =="__main__":
    p = Pool(32)
    results = p.map(TreeAve,IDquench)
    p.close()
    p.join()

    
for i in np.arange(len(results)):
    for field in Fields:
        res[field] = np.append(res[field],results[i][field])
        
np.savez('/srv/cosmdatb/erebos/yingzhong/selData/data/History_quench.npz',
         M_gas_max = res['M_gas_max'],Age_gas_max = res['Age_gas_max'],M_star_max = res['M_star_max'],
         M_bh_max = res['M_bh_max'],M_gas_befor = res['M_gas_befor'],M_star_befor = res['M_star_befor'],
         M_bh_befor = res['M_bh_befor'],M_gas_after = res['M_gas_after'],
         M_star_after = res['M_star_after'],M_bh_after = res['M_bh_after'],
         PresentSubHaloID = res['PresentSubHaloID'],BirthTime = res['BirthTime'])
