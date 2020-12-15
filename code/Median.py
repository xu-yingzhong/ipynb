import numpy as np
import illustris_python1 as il
from scipy import integrate

'''get the median number of the StarMass, bhMass and Re of some galaxies at every snapshot
   unit : M: 10**10 M_{\odot}/h Re: ckpc/h 
'''

def Median(ID):
    
    def AgeRedshift(x):
    
        h =  0.6774
        Omega_L0 = 0.6911
        Omega_r0 = 0       #4.2 * 10**(-5) / h**2
        Omega_m0 = 0.3089

        def E(z):
            return np.sqrt(Omega_L0 - Omega_r0*(1+z)**2 + Omega_m0 * (1+z)**3 + Omega_r0*(1+z)**4)
            #return np.sqrt(Omega_L0 + Omega_m0 * (1+z)**3 + Omega_r0*(1+z)**4)
        def f(z):
            return 9.785 / (h*E(z)*(1+z))
    

        return integrate.quad(lambda z:f(z),x,np.inf)[0]

    result   = {'SnapNum':np.zeros(0,dtype=int),'MstarMedian':np.zeros(0),
               'ReMedian':np.zeros(0),'MbhMedian':np.zeros(0)}
    Names    = ['MstarMedian','MbhMedian','ReMedian']
    Arr      = {} 
    Age      = np.zeros(0)
    basePath = '/srv/cosmdatb/erebos/peng/sims.TNG/TNG300-1/output'
    snapNum  = 99
    Fields   = ['SnapNum','SubhaloMassType','SubhaloBHMass','SubhaloHalfmassRadType']
    SnapRedshift = dict(np.load('/srv/cosmdatb/erebos/yingzhong/selData/data/SnapRedshift.npz'))['Redshift']
    for name in Names:
        Arr[name] = -1*np.ones(100)
    
    '''vstack every tree'''
    for Id in ID:
        tree        = il.sublink.loadTree(basePath, snapNum, id=Id, fields=Fields, onlyMPB=True)
        TreeArr     = {}
        insertIndex = np.zeros(0,dtype=int)
        
        for i in np.arange(3):
            if i == 1:
                dataArr = tree[Fields[i+1]]
            else:
                dataArr = tree[Fields[i+1]][:,4]
                
            TreeArr[Names[i]] = dataArr     

        for i in np.arange(len(tree['SnapNum'])-1):
            if tree['SnapNum'][i]-tree['SnapNum'][i+1]==2:
                insertIndex = np.append(insertIndex,i+1)
                
        if len(insertIndex)!=0:
            for i in np.arange(3):
                if i == 1:
                    dataArr = tree[Fields[i+1]]
                else:
                    dataArr = tree[Fields[i+1]][:,4]
                    
                TreeArr[Names[i]] = np.append(dataArr[0:insertIndex[0]],-1)
                for j in np.arange(len(insertIndex)-1):
                    intermediaArr     = np.append(dataArr[insertIndex[j]:insertIndex[j+1]],-1)
                    TreeArr[Names[i]] = np.append(TreeArr[Names[i]],intermediaArr)
                TreeArr[Names[i]] = np.append(TreeArr[Names[i]],dataArr[insertIndex[-1]:])

        add = -1*np.ones(100-len(TreeArr['MstarMedian']))
        for i in np.arange(3):
            TreeArr[Names[i]] = np.append(TreeArr[Names[i]],add)     
            Arr[Names[i]]     = np.vstack((Arr[Names[i]],TreeArr[Names[i]]))
    '''get the median'''
    for i in np.arange(100):
        for name in Names:
            Poindex = np.where(Arr[name][:,i]>=0)
            RealArr = Arr[name][:,i][Poindex]
            if len(RealArr) == 0:
                break
            if name == 'MstarMedian':
                result['SnapNum']  = np.append(result['SnapNum'],99-i)
            result[name] = np.append(result[name],np.median(RealArr))
            
    '''unification'''
    result['MstarMedian'] = result['MstarMedian']/result['MstarMedian'][0]
    result['MbhMedian']   = result['MbhMedian']/result['MbhMedian'][0]
    
    ScaleFactor        = 1/(1+SnapRedshift[result['SnapNum']])
    result['ReMedian'] = result['ReMedian'] * ScaleFactor
    result['ReMedian'] = result['ReMedian']/result['ReMedian'][0]
    
    for i in np.arange(len(result['SnapNum'])):
        Age = np.append(Age,AgeRedshift(SnapRedshift[result['SnapNum']][i]))
    
    result['LookBackTime'] = AgeRedshift(0) - Age

    return result

fMmass  = dict(np.load('/srv/cosmdatb/erebos/yingzhong/selData/data/haveBH.npz'))
f       = dict(np.load('/srv/cosmdatb/erebos/yingzhong/selData/data/haveBH_025.npz'))
index   = np.intersect1d(np.where(np.array(fMmass['sSFR']) > -11),np.where(np.array(fMmass['ratioS_T']) <0.25))
ID1     = f['subHaloID']
ID2     = fMmass['subHaloID'][index]
result1 = Median(ID=ID1)
result2 = Median(ID=ID2)

np.savez('/srv/cosmdatb/erebos/yingzhong/selData/data/Median_quench.npz',
        SnapNum = result1['SnapNum'],MstarMedian = result1['MstarMedian'],ReMedian = result1['ReMedian'],
        MbhMedian = result1['MbhMedian'],LookBackTime = result1['LookBackTime'])

np.savez('/srv/cosmdatb/erebos/yingzhong/selData/data/Median_form.npz',
        SnapNum = result2['SnapNum'],MstarMedian = result2['MstarMedian'],ReMedian = result2['ReMedian'],
        MbhMedian = result2['MbhMedian'],LookBackTime = result2['LookBackTime'])

