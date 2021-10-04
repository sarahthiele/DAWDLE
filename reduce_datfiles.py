from funcs_v1 import getfiles
import tqdm
import argparse
import pandas as pd

def reduce_data(pathold, pathnew, filename, label):
    bpp = pd.read_hdf(pathold+filename, key='bpp')
    
    init = bpp.groupby('bin_num').first()
    RLOFsep = bpp.loc[bpp.evol_type==3].groupby('bin_num').first()
    CEsep = bpp.loc[bpp.evol_type==7].groupby('bin_num').first()
    
    newbpp = init.append(RLOFsep)
    newbpp = newbpp.append(CEsep)
    newbpp = newbpp.sort_values(by=['bin_num', 'tphys'])
    newbpp.to_hdf(pathnew + 'new_' + filename, key='bpp')
    
    
    if label == '11_11':
        conv = bpp.loc[(bpp.kstar_1==11)&(bpp.kstar_2==11)].groupby('bin_num').first()  
    elif label == '12':
        conv = bpp.loc[(bpp.kstar_1==12)&(bpp.kstar_2.isin([10,11,12]))].groupby('bin_num').first() 
    else:
        conv = pd.read_hdf(pathold+filename, key='conv')
    conv.to_hdf(pathnew + 'reduced_' + filename, key='conv')
    conv = []
    bpp = []
    
    try:
        mass_binaries = pd.read_hdf(pathold+filename, key='mass_stars')
    except:
        print('m_binaries key')
        mass_binaries = pd.read_hdf(pathold+filename, key='mass_binaries')
        
    mass_binaries.to_hdf(pathnew + 'reduced_' + filename, key='mass_stars')
    
    return 'reduced_' + filename

parser = argparse.ArgumentParser()
parser.add_argument("--dat_path", default="./", type=str)
parser.add_argument("--dat_path_new", default="./", type=str)

args = parser.parse_args()

kstar1_list = ['10', '11', '11', '12']
kstar2_list = ['10', '10', '11', '12_10']

for kstar1, kstar2 in tqdm.tqdm(zip(kstar1_list[2:], kstar2_list[2:])):
    fnames, label = getfiles(kstar1=kstar1, kstar2=kstar2)
    for f in tqdm.tqdm(fnames):
        newf = reduce_data(pathold=args.dat_path, pathnew=args.dat_path_new, filename=f, label=label)

