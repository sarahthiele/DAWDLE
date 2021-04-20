#===================================================================================
# Old metallicity/binary fraction params:
#===================================================================================

met_arr = np.logspace(np.log10(1e-4), np.log10(0.03), 13)
met_arr = np.round(met_arr, 8)
met_arr = np.append(0.0, met_arr)

binfracs = np.array([0.4847, 0.4713, 0.458 , 0.4446, 0.4312, 0.4178, 
                     0.4044, 0.3717, 0.3309, 0.2901, 0.2493, 0.2085, 0.1677])

ratios = np.array([0.6786, 0.7171, 0.756 , 0.7971, 0.8427, 0.8904, 0.9396, 
                   1.0792, 1.2913, 1.5631, 1.9211, 2.4244, 3.1692])

ratio_05 = 0.6389

#===================================================================================
# Old file functions:
#===================================================================================

def getfiles_He_He(path):
    filename1 = path + 'dat_kstar1_10_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.0001.h5'
    filename2 = path + 'dat_kstar1_10_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.00016085.h5'
    filename3 = path + 'dat_kstar1_10_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.00025873.h5'
    filename4 = path + 'dat_kstar1_10_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.00041618.h5'
    filename5 = path + 'dat_kstar1_10_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.00066943.h5'
    filename6 = path + 'dat_kstar1_10_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.0010768.h5'
    filename7 = path + 'dat_kstar1_10_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.00173205.h5'
    filename8 = path + 'dat_kstar1_10_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.00278604.h5'
    filename9 = path + 'dat_kstar1_10_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.0044814.h5'
    filename10 = path + 'dat_kstar1_10_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.00720843.h5'
    filename11 = path + 'dat_kstar1_10_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.01159492.h5'
    filename12 = path + 'dat_kstar1_10_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.01865067.h5'
    filename13 = path + 'dat_kstar1_10_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.03.h5'

    filename_list = [filename1, filename2, filename3, filename4, filename5,
                     filename6, filename7, filename8, filename9, filename10,
                     filename11, filename12, filename13]
    label = '10_10'
    return filename_list, label

def getfiles_CO_He(path):
    filename1 = path + 'dat_kstar1_11_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.0001.h5'
    filename2 = path + 'dat_kstar1_11_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.00016085.h5'
    filename3 = path + 'dat_kstar1_11_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.00025873.h5'
    filename4 = path + 'dat_kstar1_11_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.00041618.h5'
    filename5 = path + 'dat_kstar1_11_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.00066943.h5'
    filename6 = path + 'dat_kstar1_11_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.0010768.h5'
    filename7 = path + 'dat_kstar1_11_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.00173205.h5'
    filename8 = path + 'dat_kstar1_11_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.00278604.h5'
    filename9 = path + 'dat_kstar1_11_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.0044814.h5'
    filename10 = path + 'dat_kstar1_11_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.00720843.h5'
    filename11 = path + 'dat_kstar1_11_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.01159492.h5'
    filename12 = path + 'dat_kstar1_11_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.01865067.h5'
    filename13 = path + 'dat_kstar1_11_kstar2_10_SFstart_13700.0_SFduration_0.0_metallicity_0.03.h5'
    
    filename_list = [filename1, filename2, filename3, filename4, filename5,
                     filename6, filename7, filename8, filename9, filename10,
                     filename11, filename12, filename13]
    label = '11_10'
    return filename_list, label

def getfiles_CO_CO(path):
    filename1 = path + 'dat_kstar1_11_kstar2_11_SFstart_13700.0_SFduration_0.0_metallicity_0.0001.h5'
    filename2 = path + 'dat_kstar1_11_kstar2_11_SFstart_13700.0_SFduration_0.0_metallicity_0.00016085.h5'
    filename3 = path + 'dat_kstar1_11_kstar2_11_SFstart_13700.0_SFduration_0.0_metallicity_0.00025873.h5'
    filename4 = path + 'dat_kstar1_11_kstar2_11_SFstart_13700.0_SFduration_0.0_metallicity_0.00041618.h5'
    filename5 = path + 'dat_kstar1_11_kstar2_11_SFstart_13700.0_SFduration_0.0_metallicity_0.00066943.h5'
    filename6 = path + 'dat_kstar1_11_kstar2_11_SFstart_13700.0_SFduration_0.0_metallicity_0.0010768.h5'
    filename7 = path + 'dat_kstar1_11_kstar2_11_SFstart_13700.0_SFduration_0.0_metallicity_0.00173205.h5'
    filename8 = path + 'dat_kstar1_11_kstar2_11_SFstart_13700.0_SFduration_0.0_metallicity_0.00278604.h5'
    filename9 = path + 'dat_kstar1_11_kstar2_11_SFstart_13700.0_SFduration_0.0_metallicity_0.0044814.h5'
    filename10 = path + 'dat_kstar1_11_kstar2_11_SFstart_13700.0_SFduration_0.0_metallicity_0.00720843.h5'
    filename11 = path + 'dat_kstar1_11_kstar2_11_SFstart_13700.0_SFduration_0.0_metallicity_0.01159492.h5'
    filename12 = path + 'dat_kstar1_11_kstar2_11_SFstart_13700.0_SFduration_0.0_metallicity_0.01865067.h5'
    filename13 = path + 'dat_kstar1_11_kstar2_11_SFstart_13700.0_SFduration_0.0_metallicity_0.03.h5'
    
    filename_list = [filename1, filename2, filename3, filename4, filename5,
                     filename6, filename7, filename8, filename9, filename10,
                     filename11, filename12, filename13]
    label = '11_11'
    return filename_list, label

def getfiles_ONe(path):
    filename1 = path + 'dat_kstar1_12_kstar2_10_12_SFstart_13700.0_SFduration_0.0_metallicity_0.0001.h5'
    filename2 = path + 'dat_kstar1_12_kstar2_10_12_SFstart_13700.0_SFduration_0.0_metallicity_0.00016085.h5'
    filename3 = path + 'dat_kstar1_12_kstar2_10_12_SFstart_13700.0_SFduration_0.0_metallicity_0.00025873.h5'
    filename4 = path + 'dat_kstar1_12_kstar2_10_12_SFstart_13700.0_SFduration_0.0_metallicity_0.00041618.h5'
    filename5 = path + 'dat_kstar1_12_kstar2_10_12_SFstart_13700.0_SFduration_0.0_metallicity_0.00066943.h5'
    filename6 = path + 'dat_kstar1_12_kstar2_10_12_SFstart_13700.0_SFduration_0.0_metallicity_0.0010768.h5'
    filename7 = path + 'dat_kstar1_12_kstar2_10_12_SFstart_13700.0_SFduration_0.0_metallicity_0.00173205.h5'
    filename8 = path + 'dat_kstar1_12_kstar2_10_12_SFstart_13700.0_SFduration_0.0_metallicity_0.00278604.h5'
    filename9 = path + 'dat_kstar1_12_kstar2_10_12_SFstart_13700.0_SFduration_0.0_metallicity_0.0044814.h5'
    filename10 = path + 'dat_kstar1_12_kstar2_10_12_SFstart_13700.0_SFduration_0.0_metallicity_0.00720843.h5'
    filename11 = path + 'dat_kstar1_12_kstar2_10_12_SFstart_13700.0_SFduration_0.0_metallicity_0.01159492.h5'
    filename12 = path + 'dat_kstar1_12_kstar2_10_12_SFstart_13700.0_SFduration_0.0_metallicity_0.01865067.h5'
    filename13 = path + 'dat_kstar1_12_kstar2_10_12_SFstart_13700.0_SFduration_0.0_metallicity_0.03.h5'

    filename_list = [filename1, filename2, filename3, filename4, filename5,
                     filename6, filename7, filename8, filename9, filename10,
                     filename11, filename12, filename13]
    label = '12'
    return filename_list, label

def galaxy_files_10_10_var():
    files = ['final_galaxy_10_10_0.0001_0.4847.hdf',
             'final_galaxy_10_10_0.00016085_0.4713.hdf',
             'final_galaxy_10_10_0.00025873_0.458.hdf',
             'final_galaxy_10_10_0.00041618_0.4446.hdf',
             'final_galaxy_10_10_0.00066943_0.4312.hdf',
             'final_galaxy_10_10_0.0010768_0.4178.hdf',
             'final_galaxy_10_10_0.00173205_0.4044.hdf',
             'final_galaxy_10_10_0.00278604_0.3717.hdf',
             'final_galaxy_10_10_0.0044814_0.3309.hdf',
             'final_galaxy_10_10_0.00720843_0.2901.hdf',
             'final_galaxy_10_10_0.01159492_0.2493.hdf',
             'final_galaxy_10_10_0.01865067_0.2085.hdf',
             'final_galaxy_10_10_0.03_0.1677.hdf']
    return files

def galaxy_files_11_10_var():
    files = ['final_galaxy_11_10_0.0001_0.4847.hdf',
             'final_galaxy_11_10_0.00016085_0.4713.hdf',
             'final_galaxy_11_10_0.00025873_0.458.hdf',
             'final_galaxy_11_10_0.00041618_0.4446.hdf',
             'final_galaxy_11_10_0.00066943_0.4312.hdf',
             'final_galaxy_11_10_0.0010768_0.4178.hdf',
             'final_galaxy_11_10_0.00173205_0.4044.hdf',
             'final_galaxy_11_10_0.00278604_0.3717.hdf',
             'final_galaxy_11_10_0.0044814_0.3309.hdf',
             'final_galaxy_11_10_0.00720843_0.2901.hdf',
             'final_galaxy_11_10_0.01159492_0.2493.hdf',
             'final_galaxy_11_10_0.01865067_0.2085.hdf',
             'final_galaxy_11_10_0.03_0.1677.hdf']
    return files

def galaxy_files_11_11_var():
    files = ['final_galaxy_11_11_0.0001_0.4847.hdf',
             'final_galaxy_11_11_0.00016085_0.4713.hdf',
             'final_galaxy_11_11_0.00025873_0.458.hdf',
             'final_galaxy_11_11_0.00041618_0.4446.hdf',
             'final_galaxy_11_11_0.00066943_0.4312.hdf',
             'final_galaxy_11_11_0.0010768_0.4178.hdf',
             'final_galaxy_11_11_0.00173205_0.4044.hdf',
             'final_galaxy_11_11_0.00278604_0.3717.hdf',
             'final_galaxy_11_11_0.0044814_0.3309.hdf',
             'final_galaxy_11_11_0.00720843_0.2901.hdf',
             'final_galaxy_11_11_0.01159492_0.2493.hdf',
             'final_galaxy_11_11_0.01865067_0.2085.hdf',
             'final_galaxy_11_11_0.03_0.1677.hdf']
    return files

def galaxy_files_12_var():
    files = ['final_galaxy_12_0.0001_0.4847.hdf',
             'final_galaxy_12_0.00016085_0.4713.hdf',
             'final_galaxy_12_0.00025873_0.458.hdf',
             'final_galaxy_12_0.00041618_0.4446.hdf',
             'final_galaxy_12_0.00066943_0.4312.hdf',
             'final_galaxy_12_0.0010768_0.4178.hdf',
             'final_galaxy_12_0.00173205_0.4044.hdf',
             'final_galaxy_12_0.00278604_0.3717.hdf',
             'final_galaxy_12_0.0044814_0.3309.hdf',
             'final_galaxy_12_0.00720843_0.2901.hdf',
             'final_galaxy_12_0.01159492_0.2493.hdf',
             'final_galaxy_12_0.01865067_0.2085.hdf',
             'final_galaxy_12_0.03_0.1677.hdf']
    return files

def galaxy_files_10_10_05():
    files = ['final_galaxy_10_10_0.0001_0.5.hdf',
             'final_galaxy_10_10_0.00016085_0.5.hdf',
             'final_galaxy_10_10_0.00025873_0.5.hdf',
             'final_galaxy_10_10_0.00041618_0.5.hdf',
             'final_galaxy_10_10_0.00066943_0.5.hdf',
             'final_galaxy_10_10_0.0010768_0.5.hdf',
             'final_galaxy_10_10_0.00173205_0.5.hdf',
             'final_galaxy_10_10_0.00278604_0.5.hdf',
             'final_galaxy_10_10_0.0044814_0.5.hdf',
             'final_galaxy_10_10_0.00720843_0.5.hdf',
             'final_galaxy_10_10_0.01159492_0.5.hdf',
             'final_galaxy_10_10_0.01865067_0.5.hdf',
             'final_galaxy_10_10_0.03_0.5.hdf']
    return files

def galaxy_files_11_10_05():
    files = ['final_galaxy_11_10_0.0001_0.5.hdf',
             'final_galaxy_11_10_0.00016085_0.5.hdf',
             'final_galaxy_11_10_0.00025873_0.5.hdf',
             'final_galaxy_11_10_0.00041618_0.5.hdf',
             'final_galaxy_11_10_0.00066943_0.5.hdf',
             'final_galaxy_11_10_0.0010768_0.5.hdf',
             'final_galaxy_11_10_0.00173205_0.5.hdf',
             'final_galaxy_11_10_0.00278604_0.5.hdf',
             'final_galaxy_11_10_0.0044814_0.5.hdf',
             'final_galaxy_11_10_0.00720843_0.5.hdf',
             'final_galaxy_11_10_0.01159492_0.5.hdf',
             'final_galaxy_11_10_0.01865067_0.5.hdf',
             'final_galaxy_11_10_0.03_0.5.hdf']
    return files

def galaxy_files_11_11_05():
    files = ['final_galaxy_11_11_0.0001_0.5.hdf',
             'final_galaxy_11_11_0.00016085_0.5.hdf',
             'final_galaxy_11_11_0.00025873_0.5.hdf',
             'final_galaxy_11_11_0.00041618_0.5.hdf',
             'final_galaxy_11_11_0.00066943_0.5.hdf',
             'final_galaxy_11_11_0.0010768_0.5.hdf',
             'final_galaxy_11_11_0.00173205_0.5.hdf',
             'final_galaxy_11_11_0.00278604_0.5.hdf',
             'final_galaxy_11_11_0.0044814_0.5.hdf',
             'final_galaxy_11_11_0.00720843_0.5.hdf',
             'final_galaxy_11_11_0.01159492_0.5.hdf',
             'final_galaxy_11_11_0.01865067_0.5.hdf',
             'final_galaxy_11_11_0.03_0.5.hdf']
    return files

def galaxy_files_12_05():
    files = ['final_galaxy_12_0.0001_0.5.hdf',
             'final_galaxy_12_0.00016085_0.5.hdf',
             'final_galaxy_12_0.00025873_0.5.hdf',
             'final_galaxy_12_0.00041618_0.5.hdf',
             'final_galaxy_12_0.00066943_0.5.hdf',
             'final_galaxy_12_0.0010768_0.5.hdf',
             'final_galaxy_12_0.00173205_0.5.hdf',
             'final_galaxy_12_0.00278604_0.5.hdf',
             'final_galaxy_12_0.0044814_0.5.hdf',
             'final_galaxy_12_0.00720843_0.5.hdf',
             'final_galaxy_12_0.01159492_0.5.hdf',
             'final_galaxy_12_0.01865067_0.5.hdf',
             'final_galaxy_12_0.03_0.5.hdf']
    return files

#===================================================================================
# Old LISA & GW functions:
#===================================================================================

def chirpmass(pop):
    '''
    Calculates chirp mass for circular binaries with GW frequency f_gw = 2f_orb,
    from GW lecture paper.
    
    Returns mass in kilograms.
    '''
    m1 = pop.mass_1
    m2 = pop.mass_2
    M = (m1 * m2) ** (3 / 5) / (m1 + m2) ** (1 / 5)
    return M * M_sol

def lisa_PSD():
    from scipy.interpolate import interp1d
    '''Computes LISA power spectral density according to `
    Cornish and Robson 2018 <https://arxiv.org/pdf/1803.01944.pdf>'
    without the Galactic foreground
    
    Parameters
    ----------
    none
    Returns
    -------
    LISA_hc : interpolation of LISA sensitivity curve
    '''
    freq = np.logspace(-9,1,10000)
    # note: freq [Hz], L_arm [m], S_n [Hz^-0.5]
    L_arm = 2.5e9
    f_star = 19.09*1e-3
    P_oms = (1.5e-11)**2*(1. + (2.0e-3/freq)**4)
    P_acc = (3.0e-15)**2*(1. + (0.4e-3/freq)**2)*(1. + (freq/(8.0e-3))**4)
    P_n = (P_oms + 2.*(1. + np.cos(freq/f_star)**2)*P_acc/(2.*np.pi*freq)**4)/L_arm**2
    R = 3./10./(1. + 6./10.*(freq/f_star)**2)
    S_n = (P_n/R)
    LISA_psd = interp1d(freq, S_n)
    return LISA_psd

LISA_psd = lisa_PSD()

def h_2(pop):
    '''
    Calculates the dimensionless GW strain of the 2nd harmonic.
    
    Formulas taken from Breivik(2019).
    '''
    n = 2
    Mc = chirpmass(pop)
    porb = pop.porb_f * 24 * 3600
    f_orb = 1 / porb
    DL = pop.dist_sun * m_kpc
    hns_1 = 128 / 5 * (G * Mc) ** (10 / 3) / c ** 8
    hns_2 = (np.pi * f_orb) ** (4/3) / DL ** 2 / n **2
    hns = hns_1 * hns_2
    hn = np.sqrt(hns)
    return hn

def ASD_2(pop):
    '''
    Calculates the amplitude spectral density for a stationary
    source. It assumes a observation time of 4 years.
    
    Formula taken from Breivik (2019).
    '''
    h2 = h_2(pop)
    T_obs = 4 * 365 * 24 * 3600
    ASD = h2 * np.sqrt(T_obs)
    return ASD

def SNR_circ(pop):
    '''
    Calculates the SNR for a stationary, circular
    source. Formula taken from Breivik(2019).
    '''
    LISA_psd = lisa_PSD()
    ASD = ASD_2(pop)
    PSD = LISA_psd(pop.f_gw)
    SNRs = ASD ** 2 / PSD
    return np.sqrt(SNRs)

def chirp_circ(pop):
    '''
    Calculates the chirp of a binary system for a given n harmonic. For circular binaries 
    with GW frequency f_gw = 2f_orb, F(e) is unity and the chirp only depends on
    mass and final orbital frequency.
    '''
    M = chirpmass(pop)
    f = pop.f_gw
    f_orb = 1 / (pop.porb_f * 24 * 3600)
    f_dot = 96 / (5 * np.pi) * (G * M)**(5/3) / c**5 * (2*np.pi*f_orb)**(11/3)
    return f_dot

