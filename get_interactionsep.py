import pandas as pd
import numpy as np
from funcs_v1 import getfiles
import tqdm

met_arr = np.logspace(np.log10(1e-4), np.log10(0.03), 15)
met_arr = np.round(met_arr, 8)
met_arr = np.append(0.0, met_arr)

binfracs = np.array([0.4847, 0.4732, 0.4618, 0.4503, 0.4388, 
                     0.4274, 0.4159, 0.4044, 0.3776, 0.3426, 
                     0.3076, 0.2726, 0.2376, 0.2027, 0.1677])

ratios = np.array([0.68, 0.71, 0.74, 0.78, 0.82, 
                   0.86, 0.9, 0.94, 1.05, 1.22, 
                   1.44, 1.7 , 2.05, 2.51, 3.17])

ratio_05 = 0.64


def get_intersep(pathtodat, pathtoLband, verbose=True):
    def intersep(pathtodat, datfile, pathtoLband, toname, i, label, binfrac, verbose=verbose):
        data = pd.DataFrame() 
        columns=['bin_num', 'FIRE_index', 'met', 'rad_1', 'rad_2', 'CEsep', 'CEtime', 'RLOFsep', 'RLOFtime'] 
        data.to_hdf(toname, key='data', format='t', append=True) 

        if verbose:
            print('\n{}'.format(i))

        Z = met_arr[i+1] 

        if verbose:
            print('Z: ', Z) 
            print('binfrac: ', binfrac) 

        Lbandfile = pathtoLband + 'Lband_{}_{}_{}.hdf'.format(label, Z, binfrac)
        if verbose:
            print('Lbandfile: ' + Lbandfile)
        Lband = pd.read_hdf(Lbandfile, key='Lband').sort_values('bin_num') 
        data = Lband[['bin_num', 'FIRE_index', 'met', 'rad_1', 'rad_2']] 

        if verbose:
            print('dat file: ' + datfile)
        dat = pd.read_hdf(pathtodat+datfile, key='bpp') 
        dat = dat[['tphys', 'evol_type', 'sep', 'bin_num']]

        RLOFsep = dat.loc[dat.evol_type==3].groupby('bin_num', as_index=False).first()
        RLOFsep = RLOFsep.loc[RLOFsep.bin_num.isin(data.bin_num)]
        data_RLOF = data.loc[data.bin_num.isin(RLOFsep.bin_num)]
        RLOFsep['weights'] = data_RLOF.bin_num.value_counts().sort_index().values
        
        CEsep = dat.loc[dat.evol_type==7].groupby('bin_num', as_index=False).first()
        CEsep = CEsep.loc[CEsep.bin_num.isin(data.bin_num)]
        data_CE = data.loc[data.bin_num.isin(CEsep.bin_num)]
        CEsep['weights'] = data_CE.bin_num.value_counts().sort_index().values

        dat = []
        data_RLOF = []
        data_CE = []

        data['CEsep'] = np.repeat(CEsep['sep'], CEsep['weights']).values
        data['CEtime'] = np.repeat(CEsep['tphys'], CEsep['weights']).values
        data['RLOFsep'] = np.repeat(RLOFsep['sep'], RLOFsep['weights']).values
        data['RLOFtime'] = np.repeat(RLOFsep['tphys'], RLOFsep['weights']).values

        Ntot = len(data) 
        
        if verbose:
            print('Ntot: ', Ntot) 
        N = 0 
        j = 0 
        jlast = int(1e5) 
        while j < Ntot: 
            if verbose:
                print('j: ', j) 
                print('jlast: ', jlast) 
            data[j:jlast].to_hdf(toname, key='data', format='t', append=True) 
            N += len(data[j:jlast]) 
            j += 1e5 
            j = int(j) 
            jlast += 1e5 
            if jlast > Ntot: 
                jlast = Ntot 
            jlast = int(jlast) 
        if verbose:
            if N != Ntot: 
                print('loop is wrong') 
            else: 
                print('wrote to hdf successfully') 

        return


    # FZ:

    kstar1_list = ['10', '11', '11', '12']
    kstar2_list = ['10', '10', '11', '10_12']
    for kstar1, kstar2 in zip(kstar1_list, kstar2_list):
        files, label = getfiles(kstar1, kstar2)
        for f, i in tqdm.tqdm(zip(files, range(len(files))), total=len(files)):
            if verbose:
                print('i = {}'.format(i))
            binfrac = binfracs[i]
            met = met_arr[i+1]
            toname = 'intersepfiles/{}_intersep_FZ.hdf'.format(label)
            intersep(pathtodat, f, pathtoLband, toname, i, label, binfrac, verbose)
            
            if verbose:
                print('i = {}'.format(i))
            binfrac = 0.5
            met = met_arr[i+1]
            toname = 'intersepfiles/{}_intersep_F50.hdf'.format(label)
            intersep(pathtodat, f, pathtoLband, toname, i, label, binfrac, verbose)
        
    return