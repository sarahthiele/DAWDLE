from funcs_v1 import *

pathold = '../datfiles/'
pathnew = 'reduced_datfiles/'

def reduce_data(pathold, pathnew, filename, label):
    print(filename)
    bpp = pd.read_hdf(pathold+filename, key='bpp')
    
    init = bpp.groupby('bin_num').first()
    #form1 = bpp.loc[bpp.kstar_1.isin([10, 11, 12])].groupby('bin_num').first()
    #form2 = bpp.loc[bpp.kstar_2.isin([10, 11, 12])].groupby('bin_num').first()
    RLOFsep = bpp.loc[bpp.evol_type==3].groupby('bin_num').first()
    CEsep = bpp.loc[bpp.evol_type==7].groupby('bin_num').first()
    
    #newbpp = init.append(form1)
    #newbpp = newbpp.append(form2)
    newbpp = newbpp.append(RLOFsep)
    newbpp = newbpp.append(CEsep)
    newbpp = newbpp.sort_values(by=['bin_num', 'tphys'])
    newbpp.to_hdf(pathnew + 'new_' + filename, key='bpp')
    
    if label == '10_10' or label == '11_10':
        conv = pd.read_hdf(pathold+filename, key='conv')
    elif label == '11_11':
        conv = bpp.loc[(bpp.kstar_1==11)&(bpp.kstar_2==11)].groupby('bin_num').first()  
    elif label == '12':
         conv = bpp.loc[(bpp.kstar_1==12)&(bpp.kstar_2.isin([10,11,12]))].groupby('bin_num').first() 
    conv.to_hdf(pathnew + 'reduced_' + filename, key='conv')
            
    try:
        mass_binaries = pd.read_hdf(pathold+filename, key='mass_stars')
    except:
        print('m_binaries key')
        mass_binaries = pd.read_hdf(pathold+filename, key='mass_binaries')
        
    mass_binaries.to_hdf(pathnew + 'reduced_' + filename, key='mass_stars')
    
    return 'reduced_' + filename



fnames, label = getfiles_He_He('')
for f in fnames:
    newf = reduce_data(pathold, pathnew, f, label)

fnames, label = getfiles_CO_He('')
for f in fnames:
    newf = reduce_data(pathold, pathnew, f, label)

fnames, label = getfiles_CO_CO('')
for f in fnames:
    newf = reduce_data(pathold, pathnew, f, label)

fnames, label = getfiles_ONe('')
for f in fnames:
    newf = reduce_data(pathold, pathnew, f, label)
