from funcs_v1 import Lband_files, galaxy_files
import numpy as np
import pandas as pd

def get_numLISA(pathtoLband, Lbandfile, FIREmin=0.00015, FIREmax=13.346, Z_sun=0.02):
    num = 30
    met_bins = np.logspace(np.log10(FIREmin), np.log10(FIREmax), num)*Z_sun
    
    if Lbandfile == 'old':
        He = pd.DataFrame()
        for f in galaxy_files(kstar1='10', kstar2='10', var=True):
            He = He.append(pd.read_hdf(pathtoLband + f, key='Lband'))  
        print('finished He + He')
        COHe = pd.DataFrame()
        for f in galaxy_files(kstar1='11', kstar2='10', var=True):
            COHe = COHe.append(pd.read_hdf(pathtoLband + f, key='Lband'))
        print('finished CO + He')
        CO = pd.DataFrame()
        for f in galaxy_files(kstar1='11', kstar2='11', var=True):
            CO = CO.append(pd.read_hdf(pathtoLband + f, key='Lband'))
        print('finished CO + CO')
        ONe = pd.DataFrame()
        for f in galaxy_files(kstar1='12', kstar2='10', var=True):
            ONe = ONe.append(pd.read_hdf(pathtoLband + f, key='Lband'))
        print('finished ONe + X')
    
    elif Lbandfile == 'new':
        He = pd.DataFrame()
        for f in Lband_files(kstar1='10', kstar2='10', var=True):
            He = He.append(pd.read_hdf(pathtoLband + f, key='Lband'))
        print('finished He + He')
        COHe = pd.DataFrame()
        for f in Lband_files(kstar1='11', kstar2='10', var=True):
            COHe = COHe.append(pd.read_hdf(pathtoLband + f, key='Lband'))
        print('finished CO + He')
        CO = pd.DataFrame()
        for f in Lband_files(kstar1='11', kstar2='11', var=True):
            CO = CO.append(pd.read_hdf(pathtoLband + f, key='Lband'))
        print('finished CO + CO')
        ONe = pd.DataFrame()
        for f in Lband_files(kstar1='12', kstar2='10', var=True):
            ONe = ONe.append(pd.read_hdf(pathtoLband + f, key='Lband'))
        print('finished ONe + X')
        
    Henums, bins = np.histogram(He.met*Z_sun, bins=met_bins)
    COHenums, bins = np.histogram(COHe.met*Z_sun, bins=met_bins)
    COnums, bins = np.histogram(CO.met*Z_sun, bins=met_bins)
    ONenums, bins = np.histogram(ONe.met*Z_sun, bins=met_bins)

    numLISA_30bins = pd.DataFrame(np.array([Henums, COHenums, COnums, ONenums]).T, 
                                     columns=['He', 'COHe', 'CO', 'ONe'])

    numLISA_30bins.to_hdf('numLISA/numLISA_30bins.hdf', key='data')

    # F50:
    
    if Lbandfile == 'old':
        He05 = pd.DataFrame()
        for f in galaxy_files(kstar1='10', kstar2='10', var=False):
            He05 = He05.append(pd.read_hdf(pathtoLband + f, key='Lband'))
        print('finished He + He, F50')
        COHe05 = pd.DataFrame()
        for f in galaxy_files(kstar1='11', kstar2='10', var=False):
            COHe05 = COHe05.append(pd.read_hdf(pathtoLband + f, key='Lband'))
        print('finished CO + He, F50')
        CO05 = pd.DataFrame()
        for f in galaxy_files(kstar1='11', kstar2='11', var=False):
            CO05 = CO05.append(pd.read_hdf(pathtoLband + f, key='Lband'))
        print('finished CO + CO, F50')
        ONe05 = pd.DataFrame()
        for f in galaxy_files(kstar1='12', kstar2='10', var=False):
            ONe05 = ONe05.append(pd.read_hdf(pathtoLband + f, key='Lband'))
        print('finished ONe + X, F50') 

    elif Lbandfile == 'new':
        He05 = pd.DataFrame()
        for f in Lband_files(kstar1='10', kstar2='10', var=False):
            He05 = He05.append(pd.read_hdf(pathtoLband + f, key='Lband'))
        print('finished He + He, F50')
        COHe05 = pd.DataFrame()
        for f in Lband_files(kstar1='11', kstar2='10', var=False):
            COHe05 = COHe05.append(pd.read_hdf(pathtoLband + f, key='Lband'))
        print('finished CO + He, F50')
        CO05 = pd.DataFrame()
        for f in Lband_files(kstar1='11', kstar2='11', var=False):
            CO05 = CO05.append(pd.read_hdf(pathtoLband + f, key='Lband'))
        print('finished CO + CO, F50')
        ONe05 = pd.DataFrame()
        for f in Lband_files(kstar1='12', kstar2='10', var=False):
            ONe05 = ONe05.append(pd.read_hdf(pathtoLband + f, key='Lband'))
        print('finished ONe + X, F50')
    
    Henums05, bins = np.histogram(He05.met*Z_sun, bins=met_bins)
    COHenums05, bins = np.histogram(COHe05.met*Z_sun, bins=met_bins)
    COnums05, bins = np.histogram(CO05.met*Z_sun, bins=met_bins)
    ONenums05, bins = np.histogram(ONe05.met*Z_sun, bins=met_bins)

    numLISA_30bins_05 = pd.DataFrame(np.array([Henums05, COHenums05, COnums05, ONenums05]).T, 
                                     columns=['He', 'COHe', 'CO', 'ONe'])

    numLISA_30bins_05.to_hdf('numLISA/numLISA_30bins_05.hdf', key='data')
    
    return


