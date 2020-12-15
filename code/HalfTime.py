import numpy as np
import illustris_python1 as il
from multiprocessing import Pool
from functools import reduce


def HalfTime(id):
    
    def sort(Arr,Num):
        Arr   = np.absolute(Arr - Num/2)
        index = np.argsort(Arr)[0]
    
        return index
    
    basePath    = '/srv/cosmdatb/erebos/peng/sims.TNG/TNG300-1/output'
    snapNum     = 99
    Fields      = ['SnapNum','SubhaloMassType','SubfindID','SubhaloBHMass','SubhaloHalfmassRadType']
    Redshift    = dict(np.load('/srv/cosmdatb/erebos/yingzhong/selData/data/SnapRedshift.npz'))['Redshift']
    
    ScaleFactor = 1/(Redshift+1)
    result   = {}
    
    tree = il.sublink.loadTree(basePath, snapNum, id=id, fields=Fields, onlyMPB=True)
    
    for i in np.arange(len(tree['SnapNum'])):
        tree['SubhaloHalfmassRadType'][:,4][i] = ScaleFactor[tree['SnapNum'][i]] *  tree['SubhaloHalfmassRadType'][:,4][i]
    
    indexStar = sort(tree['SubhaloMassType'][:,4],tree['SubhaloMassType'][:,4][0])
    indexHalo = sort(tree['SubhaloMassType'][:,1],tree['SubhaloMassType'][:,1][0])
    indexBH1  = sort(tree['SubhaloMassType'][:,5],tree['SubhaloMassType'][:,5][0])
    indexBH   = sort(tree['SubhaloBHMass'],tree['SubhaloBHMass'][0])
    indexRe   = sort(tree['SubhaloHalfmassRadType'][:,4],tree['SubhaloHalfmassRadType'][:,4][0])
    
    result['RedshiftHalfStar'] = Redshift[tree['SnapNum'][indexStar]]
    result['RedshiftHalfHalo'] = Redshift[tree['SnapNum'][indexHalo]]
    result['RedshiftHalfBH']   = Redshift[tree['SnapNum'][indexBH]]
    result['RedshiftHalfBH1']  = Redshift[tree['SnapNum'][indexBH1]]
    result['RedshiftHalRe']    = Redshift[tree['SnapNum'][indexRe]]
    
    return result

fMmass = dict(np.load('/srv/cosmdatb/erebos/yingzhong/selData/data/haveBH.npz'))
f1     = dict(np.load('/srv/cosmdatb/erebos/yingzhong/selData/data/haveBH_025.npz'))

index  = reduce(np.intersect1d,
                (np.where(np.array(fMmass['sSFR']) > -11),
                 np.where(np.array(fMmass['ratioS_T']) <= 0.25)))

Indexform   = fMmass['subHaloID'][index]
Indexquench = f1['subHaloID']

form  = {'RedshiftHalfStar':np.zeros(0),
         'RedshiftHalfHalo':np.zeros(0),
         'RedshiftHalfBH':np.zeros(0),
         'RedshiftHalfBH1':np.zeros(0),
          'RedshiftHalRe':np.zeros(0)}

#quench = {'RedshiftHalfStar':np.zeros(0),
#         'RedshiftHalfHalo':np.zeros(0),
#         'RedshiftHalfBH':np.zeros(0),
#         'RedshiftHalfBH1':np.zeros(0),
#          'RedshiftHalRe':np.zeros(0)}

if __name__ =="__main__":
    p = Pool(32)
    results = p.map(HalfTime,Indexform)
    p.close()
    p.join()
    
for i in np.arange(len(results)):
    form['RedshiftHalfStar'] = np.append(form['RedshiftHalfStar'],results[i]['RedshiftHalfStar'])
    form['RedshiftHalfHalo'] = np.append(form['RedshiftHalfHalo'],results[i]['RedshiftHalfHalo'])
    form['RedshiftHalfBH']   = np.append(form['RedshiftHalfBH'],results[i]['RedshiftHalfBH'])
    form['RedshiftHalfBH1']  = np.append(form['RedshiftHalfBH1'],results[i]['RedshiftHalfBH1'])
    form['RedshiftHalRe']    = np.append(form['RedshiftHalRe'],results[i]['RedshiftHalRe'])
    
np.savez('/srv/cosmdatb/erebos/yingzhong/selData/data/Halftime_form.npz',
         RedshiftHalfStar = form['RedshiftHalfStar'],RedshiftHalfHalo=form['RedshiftHalfHalo'],
         RedshiftHalfBH = form['RedshiftHalfBH'],RedshiftHalfBH1=form['RedshiftHalfBH1'],
         RedshiftHalRe = form['RedshiftHalRe'])   
    