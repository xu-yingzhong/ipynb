import numpy as np
import illustris_python1 as il
from matplotlib import pyplot as plt
import seaborn as sns
import os
import h5py
'''==========================================================================================='''
'''
    periodic distance 
    All the  units of distance here is ckpc/h. assume the radius of the galaxy is < L/2
'''
'''==========================================================================================='''
def periodDistance(distance,L = 205000):
    if distance > L/2 :
        return  distance - L
    else:
        if distance < -L/2 :
            return distance + L
        else:
            return distance

'''==========================================================================================='''
'''
    Rodrigues' rotation matrix 
    (x,y,z), unit vector, represents the fixed axis, theta is the angle (arc)
'''
'''==========================================================================================='''
def rotationMat(x,y,z,theta):
    w = np.array([[0,-z,y],[z,0,-x],[-y,x,0]])
    I = np.identity(3)
    R = I + w*np.sin(theta) + np.dot(w,w)*(1-np.cos(theta))

    return R

'''============================================================================================'''
def Null():
    result = {}
    rfields = ['Coordinates','Velocities','cirPar','Masses','radius']
    for field in rfields:
        result[field] = np.array([])

    return result

'''============================================================================================='''
'''
    get the circularity parameter 

    ID: subhaloID
    a: 1/(1+z) scalfactor
    H_z_0: H_{z}/H_{0} the ratio of Hubble parameter
'''
'''==============================================================================================='''
def cirparameter(basePath,snapNum,ID,a,H_z_0,q=1):
    
    result1 = {}
    fields = ['Coordinates','Velocities','GFM_StellarFormationTime','Masses']

    '''load position, and vercility of stars in the subhalo'''

    parData = il.snapshot.loadSubhalo(basePath,snapNum,ID,'star',fields=fields)

    '''find the the position and velocity of the subhalo'''

    subhaloInf = il.groupcat.loadSingle(basePath,snapNum,subhaloID=ID)
    subhaloPos = subhaloInf['SubhaloPos']
    subhaloVel = subhaloInf['SubhaloVel']
    starNum = subhaloInf['SubhaloLenType'][4]
    R_M = subhaloInf['SubhaloHalfmassRadType'][4] * a

    if starNum == 0:
        result1 = Null()

        return result1
    
    '''reset the positions and verlosities '''

    for i in np.arange(starNum):
        for j in np.arange(3):
            d = parData['Coordinates'][i,j]-subhaloPos[j]  # unit ckpc/h
            cDistance = periodDistance(distance=d)         # period box
            parData['Coordinates'][i,j] = cDistance * a    # unit kpc/h
            parData['Velocities'][i,j] = parData['Velocities'][i,j] * np.sqrt(a) - subhaloVel[j]
            parData['Velocities'][i,j] += H_z_0 * parData['Coordinates'][i,j] * 0.1  # change unit to km/s

    '''slect the real star particles and radius <3*R_M'''
    obj = []
    starRad = np.zeros(starNum)
    
    for i in np.arange(starNum): #get the radii of every star and wind particle
        starRad[i] = np.linalg.norm(parData['Coordinates'][i,:])

    parData['radius'] = starRad
    
    for i in np.arange(starNum):
        if parData['GFM_StellarFormationTime'][i] <= 0 or starRad[i] >= 3*R_M :   
            obj.append(i)

    fields = ['Coordinates','Velocities','GFM_StellarFormationTime','Masses','radius']
    for field in fields:
        parData[field] = np.delete(parData[field],obj,axis = 0)

    starNum = len(parData['GFM_StellarFormationTime']) # the number of real stars within 3*R_M
    if starNum == 0:
        result1 = Null()

        return result1

    #change the unit of mass to $M_{\odot}/h$
    parData['Masses'] = parData['Masses']*10**10
    
    '''get the radius R within which the number of stars is q% of the total '''
    #radLarToSma = np.array(sorted(parData['radius'],reverse=True))
    #R_qIndex = round(starNum*(1-q)) - 1
    #if R_qIndex >= 0:
    #    R_q = radLarToSma[R_qIndex]
    #    innerIndex = np.where(parData['radius']<=R_q)
    #    for field in fields:
    #        parData[field] = parData[field][innerIndex]

    #    starNum = len(parData['GFM_StellarFormationTime'])

    '''calculate the specific angular momentum'''

    angularMom = np.zeros([starNum,3])
    angularNorm = np.zeros(starNum)
    angularTotalMom = np.zeros(3)
    for i in np.arange(starNum):
        angularMom[i,:] = np.cross(parData['Coordinates'][i,:],parData['Velocities'][i,:])
        angularTotalMom = angularMom[i,:] + angularTotalMom

    parData['angularMomentum'] = angularMom
    angularTotalMom = angularTotalMom / np.linalg.norm(angularTotalMom)

    ''' rotate '''
    e_z = np.array([0,0,1])
    rotAxis = np.cross(angularTotalMom,e_z)
    rotAxis = rotAxis/np.linalg.norm(rotAxis)
    cosTheta = np.dot(e_z,angularTotalMom)
    theta = np.arccos(cosTheta)
    R = rotationMat(rotAxis[0],rotAxis[1],rotAxis[2],theta = theta) #rotation matrix

    sFields = ['Coordinates','Velocities','angularMomentum']
    for field in sFields:
        for i in np.arange(starNum):
            parData[field][i,:] = np.dot(R,parData[field][i,:])

    ''' get the circularity parameter for every particle'''

    for i in np.arange(starNum):
        angularNorm[i] = np.linalg.norm(parData['angularMomentum'][i,:])

    cirPar = np.zeros(starNum)
    for i in np.arange(starNum):
        if angularNorm[i] == 0:
            continue
        cirPar[i] = parData['angularMomentum'][i,2]/angularNorm[i]

    parData['cirPar'] = cirPar
    
    rfields = ['Coordinates','Velocities','cirPar','Masses','radius']
    for field in rfields:
        result1[field] = parData[field]

    return result1

'''========================================================================================'''
'''
   get SFR
'''
'''========================================================================================'''
def getSFR(basePath,snapNum,ID,a):
    fields = ['Coordinates','StarFormationRate']
    '''load Data'''
    parData = il.snapshot.loadSubhalo(basePath,snapNum,ID,'gas',fields=fields)
    subhaloInf = il.groupcat.loadSingle(basePath,snapNum,subhaloID=ID)
    subhaloPos = subhaloInf['SubhaloPos']
    starNum = subhaloInf['SubhaloLenType'][0]
    R_M = subhaloInf['SubhaloHalfmassRadType'][4] * a

    if starNum == 0:
        SFR = 0

        return SFR

    '''reset the position'''
    for i in np.arange(starNum):
        for j in np.arange(3):
            d = parData['Coordinates'][i,j]-subhaloPos[j]  # unit ckpc/h
            cDistance = periodDistance(distance=d)
            parData['Coordinates'][i,j] = cDistance * a    # unit kpc/h

    '''select gas cell within 3*R_M'''
    starRad = np.zeros(starNum)
    for i in np.arange(starNum):
        starRad[i] = np.linalg.norm(parData['Coordinates'][i,:])

    innerGasId = np.where(starRad<=3*R_M)
    SFR = np.sum(parData['StarFormationRate'][innerGasId])

    return SFR

'''==========================================================================================='''
'''
    select subhalo 

    units:
      limStarMass : M_{\odot}
      sSfrLog     : log_{10}(yr^{-1})

    trigger=True/False : turn on/off the operation about selecting sSFR
'''
'''============================================================================================'''
def selectSubhalo(basePath,snapNum,a,H_z_0,ID,limStarMass,h=0.677):
    result = [0,0,0]
    SpMass = 0
    SFR = getSFR(basePath,snapNum,ID,a)
    Data = cirparameter(basePath,snapNum,ID,a,H_z_0)
    starNum = len(Data['Masses'])

    Data['Masses'] = Data['Masses']/h             # unit of mass: M_{\odot}
    TotalStellarMass = np.sum(Data['Masses'])
    
    if TotalStellarMass <= limStarMass:
        for i in np.arange(3):
            result[i] = []
            
        return result
    
    if SFR == 0:
        sSFR = -20                                # can choose another number
    else:
        sSFR = np.log10(SFR/TotalStellarMass)                   # unit of sSFR: yr^{-1}

    for i in np.arange(starNum):
        if Data['cirPar'][i]<0:
            SpMass = SpMass + Data['Masses'][i]*2
        elif Data['cirPar'][i]==0:
            SpMass = SpMass + Data['Masses'][i]
        else:
            continue
    ratioS_T = SpMass/TotalStellarMass

    result[0] = [ID]
    result[1] = [ratioS_T]
    result[2] = [sSFR]

    return result

'''==============================================================================================='''
'''
   accurate S-D decomposition
'''
'''==============================================================================================='''
def SDseparateCC(basePath,snapNum,subhaloID,a,H_z_0):
    
    result = dict(spheroid={},disc={})
    ''' load Data'''                                                    
    parData = cirparameter(basePath,snapNum,subhaloID,a,H_z_0)
    radius = parData['radius']
    cirPar = parData['cirPar']
    parNum = len(cirPar)
    index = np.arange(parNum)
    spherIndex = np.zeros(0,dtype=int)
    '''search and register the right spheroid stellar'''
    indexL = index[np.where(cirPar<0)]

    cirParR = cirPar[np.where(cirPar>0)]
    radiusR = radius[np.where(cirPar>0)]
    indexR = index[np.where(cirPar>0)]

    for i in indexL:
        if len(indexR) == 0:           #the subhalo,considered,must have S/T <0.5
            break
        
        rSubtrArrayR = radius[i]-radiusR
        rSubtrSqArrR = np.power(rSubtrArrayR,2)
        cSubtrArrayR = -1*cirPar[i]-cirParR
        cSubtrSqArrR = np.power(cSubtrArrayR,2)
        totalRadSqR = rSubtrSqArrR + cSubtrSqArrR
        nearestIndex = np.argmin(totalRadSqR)
        spherIndex = np.append(spherIndex,indexR[nearestIndex])
        cirParR = np.delete(cirParR,nearestIndex)
        radiusR = np.delete(radiusR,nearestIndex)
        indexR = np.delete(indexR,nearestIndex)

    for i in np.arange(parNum):         
        if cirPar[i]<=0 :
            spherIndex = np.append(spherIndex,i)

    discIndex = np.delete(np.arange(parNum),spherIndex,axis=0)
    #discCirpar = parData['cirPar'][discIndex]
    #discIndex = discIndex[np.where(discCirpar>=0.7)]
        
    
    rfields = ['Coordinates','Velocities','cirPar','Masses','radius']
    for field in rfields:
        result['spheroid'][field] = parData[field][spherIndex]

    result['spheroid']['Index'] = spherIndex

    for field in rfields:
        result['disc'][field] = parData[field][discIndex]

    result['disc']['Index'] = discIndex

    return result  
              
'''==============================================================================================='''
'''
   M_star 
'''
'''==============================================================================================='''
def M_star(basePath,snapNum,ID,a,h):
    fields = ['Coordinates','GFM_StellarFormationTime','Masses']
    
    '''load positions of stars in the subhalo'''
    
    parData = il.snapshot.loadSubhalo(basePath,snapNum,ID,'star',fields=fields)
    
    '''find the the position of the subhalo'''
    
    subhaloInf = il.groupcat.loadSingle(basePath,snapNum,subhaloID=ID)
    subhaloPos = subhaloInf['SubhaloPos']
    starNum = subhaloInf['SubhaloLenType'][4]
    R_M = subhaloInf['SubhaloHalfmassRadType'][4] * a      # unit kpc/h

    '''reset the positions'''

    for i in np.arange(starNum):
        for j in np.arange(3):
            d = parData['Coordinates'][i,j]-subhaloPos[j]  # unit ckpc/h
            cDistance = periodDistance(distance=d)         # period box
            parData['Coordinates'][i,j] = cDistance * a    # unit kpc/h

    '''slect the real star particles and radius <3*R_M'''

    starMass = 0
    starRad = np.zeros(starNum)

    for i in np.arange(starNum):                           #get the radii of every star and wind particle
        starRad[i] = np.linalg.norm(parData['Coordinates'][i,:])

    #parData['radius'] = starRad

    for i in np.arange(starNum):
        if parData['GFM_StellarFormationTime'][i] > 0 and starRad[i] < 3*R_M :
            starMass = starMass + parData['Masses'][i]

    starMass = (starMass * 10**10)/h

    return starMass
            
'''==============================================================================================='''
'''
   M_halo r<200,c
'''
'''==============================================================================================='''
def M_halo(basePath,snapNum,ID,a,h):
    address = basePath + '/snapdir_099/snap_099.0.hdf5'        # TNG '/snapdir_099/snap_099.0.hdf5'
    hf = h5py.File(address,'r')
    header = hf['Header']
    dmMass = (header.attrs['MassTable'][1] * 10**10)/h   # unit M_sun
    #dmNum_total  = header.attrs['NumPart_Total_HighWord'][1]*2**32 + header.attrs['NumPart_Total'][1]
    hf.close()
    '''find the id of the host halo'''
    subhaloInf = il.groupcat.loadSingle(basePath,snapNum,subhaloID=ID)
    haloId = subhaloInf['SubhaloGrNr']

    '''load positions of dm in the halo'''
    dmPosition = il.snapshot.loadHalo(basePath,snapNum,haloId,'dm',fields='Coordinates')

    '''find the the position of the halo'''
    haloInf = il.groupcat.loadSingle(basePath,snapNum,haloID=haloId)
    haloPos = haloInf['GroupPos']
    dmNum   = haloInf['GroupLenType'][1]
    R_200c  = haloInf['Group_R_Crit200']*a

    '''reset the positions'''
    for i in np.arange(dmNum):
        for j in np.arange(3):
            d = dmPosition[i,j]-haloPos[j]
            cDistance = periodDistance(distance=d)
            dmPosition[i,j] = cDistance*a

    '''slect dm with radius < R_200c'''
    N_dm = 0
    dmRad = np.zeros(dmNum)

    for i in np.arange(dmNum):
        dmRad[i] = np.linalg.norm(dmPosition[i,:])

    for i in np.arange(dmNum):
        if dmRad[i] < R_200c :
            N_dm +=1

    M_halo = N_dm * dmMass


    return M_halo

'''==============================================================================================='''
'''
   M_BH,central & id
'''
'''==============================================================================================='''
def BH(basePath,snapNum,ID,a,h):
    result = {}
    fields = ['Coordinates','ParticleIDs','BH_Mass']
    
    '''load positions of BHs in the subhalo'''
    
    parData = il.snapshot.loadSubhalo(basePath,snapNum,ID,'bh',fields=fields)
    
    '''find the the position of the subhalo'''
    
    subhaloInf = il.groupcat.loadSingle(basePath,snapNum,subhaloID=ID)
    subhaloPos = subhaloInf['SubhaloPos']
    bhNum = subhaloInf['SubhaloLenType'][5]

    if bhNum == 0:
        result['BH_CenMass'] = 0
        result['BH_ID'] = -1
        return result

    '''reset the positions'''

    for i in np.arange(bhNum):
        for j in np.arange(3):
            d = parData['Coordinates'][i,j]-subhaloPos[j]  # unit ckpc/h
            cDistance = periodDistance(distance=d)         # period box
            parData['Coordinates'][i,j] = cDistance * a    # unit kpc/h


    bhRad = np.zeros(bhNum)

    for i in np.arange(bhNum):                           #get the radii of every bh
        bhRad[i] = np.linalg.norm(parData['Coordinates'][i,:])

    bhIndex = np.argmin(bhRad)
    result['BH_CenMass'] = (parData['BH_Mass'][bhIndex] * 10**10)/h
    result['BH_ID'] = parData['ParticleIDs'][bhIndex]

    return result


'''==============================================================================================='''
'''
   plot HiH2 radial profiles
'''
'''==============================================================================================='''
def HiH2RadialProfile(ID,basePath,snapNum,a=1,h=0.677):
    Nr = 1
    Nc = 2
    proj  = ['_2d']
    col   = ['_neutral_H','_mol_GK11']
    sx = ['$HI+H_{2}$','$H_{2}\ $GK11']
    subhaloInf = il.groupcat.loadSingle(basePath,snapNum,subhaloID=ID)
    R_starHalf = (subhaloInf['SubhaloHalfmassRadType'][4] * a) / h
    R_gasHalf  = (subhaloInf['SubhaloHalfmassRadType'][0] * a) / h *0.1
    #haloId  = subhaloInf['SubhaloGrNr']
    #haloInf = il.groupcat.loadSingle(basePath,snapNum,haloID=haloId)
    #R_200c  = (haloInf['Group_R_Crit200'] * a) / h
    
    Data = {}
    with h5py.File('/srv/cosmdatb/erebos/peng/sims.TNG/TNG300-1/Supp/hih2_galaxy_099.hdf5','r') as f:
        if ID not in f['id_subhalo']:
            raise Exception('this subhalo is not included')
            
        indexR = np.where(np.array(f['id_subhalo']) == ID)[0][0]
        binsD   = np.array(f['profile_bins'][indexR,:]) - R_gasHalf 
        indexC_lim = np.min(np.where(binsD >= 0)) + 1 
        
        for i in np.arange(Nr):
            for j in np.arange(Nc):
                Data['profile_f'+col[j]+proj[i]] = np.array(f['profile_f'+col[j]+proj[i]][indexR,:indexC_lim])
                
            Data['profile_gas_rho'+proj[i]] = np.array(f['profile_gas_rho'+proj[i]][indexR,:indexC_lim])
                
        Data['profile_bins'] = np.array(f['profile_bins'][indexR,:indexC_lim])
        
    '''plot'''
    #fig = plt.figure()
    #gs = fig.add_gridspec(Nr, Nc, hspace=0, wspace=0)
    #axs = gs.subplots(sharex=True, sharey=True)
    fig, axs = plt.subplots(Nr, Nc,sharex=True, sharey=True)
    
    for i in np.arange(Nr):
        for j in np.arange(Nc):
            x = Data['profile_bins']
            if j == 0:
                y = np.log10(Data['profile_f'+col[j]+proj[i]] * Data['profile_gas_rho'+proj[i]] + 1)
            else:
                y = np.log10(Data['profile_f'+col[0]+proj[i]] * Data['profile_gas_rho'+proj[i]] * Data['profile_f'+col[j]+proj[i]] + 1)  
            
            axs[j].plot(x,y,label=proj[i].strip('_'))
            axs[j].set_title(sx[j],{'fontsize':'small'})
            axs[j].set(xlabel='R/kpc', ylabel='$\log_{10}(\Sigma / (M_{\odot}/kpc^2))$')
            axs[j].vlines([R_starHalf,R_gasHalf], 0, 1, transform=axs[j].get_xaxis_transform(), 
                            colors='r',linestyles=['dotted','dashed'])
            axs[j].legend(fontsize='small')
            
    for i in np.arange(Nr):
        for j in np.arange(Nc):
            axs[j].label_outer()
    fig.suptitle('ID:%d' %ID)
    
    plt.show()
    plt.savefig('/srv/cosmdatb/erebos/yingzhong/fig/HiH2Fig%d'%ID,dpi=500)