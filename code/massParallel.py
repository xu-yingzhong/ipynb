import illustris_python1 as il
import toycode1 as toy
import numpy as np
from multiprocessing import Pool


basePath = '/srv/cosmdatb/erebos/peng/sims.TNG/TNG300-1/output'
snapNum = 99
h = 0.677

def Func(ID):
    result = {}
    a = 1
    H_z_0 = 1
    result['subHaloID']  = ID
    result['M_star']     = toy.M_star(basePath,snapNum,ID,a,h)
    result['M_dm_200c']  = toy.M_halo(basePath,snapNum,ID,a,h)
    result['BH_ID']      = toy.BH(basePath,snapNum,ID,a,h)['BH_ID']
    result['BH_CenMass'] = toy.BH(basePath,snapNum,ID,a,h)['BH_CenMass']

    return result

idArr= np.zeros(0,dtype=int)
M_StarArr = np.zeros(0)
M_dmArr = np.zeros(0)
BH_idArr = np.zeros(0,dtype=int)
M_bhArr = np.zeros(0)

f = np.load('/srv/cosmdatb/erebos/yingzhong/selData/data/result_full_P.npz')
firstSubId = f['ID']


if __name__ =="__main__":
    p = Pool(32)
    results = p.map(Func,firstSubId)
    p.close()
    p.join()
       
for i in np.arange(len(results)):
    idArr        = np.append(idArr,results[i]['subHaloID'])
    M_StarArr    = np.append(M_StarArr,results[i]['M_star'])
    M_dmArr      = np.append(M_dmArr,results[i]['M_dm_200c'])
    BH_idArr     = np.append(BH_idArr,results[i]['BH_ID'])
    M_bhArr      = np.append(M_bhArr,results[i]['BH_CenMass'])
        
    
'''archive'''


np.savez('/srv/cosmdatb/erebos/yingzhong/selData/data/resultMass.npz',
         subHaloID = idArr,M_star=M_StarArr,M_dm_200c=M_dmArr,
         BH_ID=BH_idArr,BH_CenMass=M_bhArr)






      
