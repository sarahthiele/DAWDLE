from funcs_v1 import *

def get_intersep(pathtodat, pathtoLband):
    def intersep(datfile, pathtoLband, toname, i, label, binfrac):
        data = pd.DataFrame() 
        columns=['bin_num', 'FIRE_index', 'met', 'rad_1', 'rad_2', 'CEsep', 'CEtime', 'RLOFsep', 'RLOFtime'] 
        data.to_hdf(toname, key='data', format='t', append=True) 

        print('\n{}'.format(i))

        Z = met_arr[i+1] 

        print('Z: ', Z) 
        print('binfrac: ', binfrac) 

        Lbandfile = pathtoLband + '_{}_{}_{}.hdf'.format(label, Z, binfrac)
        print('Lbandfile: ' + Lbandfile)
        Lband = pd.read_hdf(Lbandfile, key='Lband').sort_values('bin_num') 
        data = Lband[['bin_num', 'FIRE_index', 'met', 'rad_1', 'rad_2']] 

        print('dat file: ' + datfile)
        dat = pd.read_hdf(datfile, key='bpp') 
        dat = dat[['tphys', 'evol_type', 'sep']]

        RLOFsep = dat.loc[dat.evol_type==3].groupby('bin_num').first()
        RLOFsep = RLOFsep.loc[Lband.bin_num.values]
        CEsep = dat.loc[dat.evol_type==7].groupby('bin_num').first()
        CEsep = CEsep.loc[Lband.bin_num.values]
        dat = []

        data['CEsep'] = CEsep.sep.values
        data['CEtime'] = CEsep.tphys.values
        data['RLOFsep'] = RLOFsep.sep.values
        data['RLOFtime'] = RLOFsep.tphys.values

        Ntot = len(data) 
        print('Ntot: ', Ntot) 
        N = 0 
        j = 0 
        jlast = int(1e5) 
        while j < Ntot: 
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
        if N != Ntot: 
            print('loop is wrong') 
        else: 
            print('wrote to hdf successfully') 

        return


    # FZ:

    files, label = getfiles_He_He(pathtodat)
    for i, f in enumerate(files):
        print('i = {}'.format(i))
        binfrac = binfracs[i]
        met = met_arr[i+1]
        toname = 'intersepfiles/He_He_intersep_FZ.hdf'.format(label, met, binfrac)
        intersep(f, pathtoLband, toname, i, label, binfrac)

    files, label = getfiles_CO_He(pathtodat)
    for i, f in enumerate(files):
        print('i = {}'.format(i))
        binfrac = binfracs[i]
        met = met_arr[i+1]
        toname = 'intersepfiles/CO_He_intersep_FZ.hdf'.format(label, met, binfrac)
        intersep(f, pathtoLband, toname, i, label, binfrac)

    files, label = getfiles_CO_CO(pathtodat)
    for i, f in enumerate(files):
        print('i = {}'.format(i))
        binfrac = binfracs[i]
        met = met_arr[i+1]
        toname = 'intersepfiles/CO_CO_intersep_FZ.hdf'.format(label, met, binfrac)
        intersep(f, pathtoLband, toname, i, label, binfrac)

    files, label = getfiles_ONe(pathtodat)
    for i, f in enumerate(files):
        print('i = {}'.format(i))
        binfrac = binfracs[i]
        met = met_arr[i+1]
        toname = 'intersepfiles/ONe_intersep_FZ.hdf'.format(label, met, binfrac)
        intersep(f, pathtoLband, toname, i, label, binfrac)

    # F50:

    files, label = getfiles_He_He(pathtodat)
    for i, f in enumerate(files):
        print('i = {}'.format(i))
        binfrac = 0.5
        met = met_arr[i+1]
        toname = 'intersepfiles/He_He_intersep_F50.hdf'.format(label, met, binfrac)
        intersep(f, pathtoLband, toname, i, label, binfrac)

    files, label = getfiles_CO_He(pathtodat)
    for i, f in enumerate(files):
        print('i = {}'.format(i))
        binfrac = 0.5
        met = met_arr[i+1]
        toname = 'intersepfiles/CO_He_intersep_F50.hdf'.format(label, met, binfrac)
        intersep(f, pathtoLband, toname, i, label, binfrac)

    files, label = getfiles_CO_CO(pathtodat)
    for i, f in enumerate(files):
        print('i = {}'.format(i))
        binfrac = 0.5
        met = met_arr[i+1]
        toname = 'intersepfiles/CO_CO_intersep_F50.hdf'.format(label, met, binfrac)
        intersep(f, pathtoLband, toname, i, label, binfrac)

    files, label = getfiles_ONe(pathtodat)
    for i, f in enumerate(files):
        print('i = {}'.format(i))
        binfrac = 0.5
        met = met_arr[i+1]
        toname = 'intersepfiles/ONe_intersep_F50.hdf'.format(label, met, binfrac)
        intersep(f, pathtoLband, toname, i, label, binfrac)
        
    return