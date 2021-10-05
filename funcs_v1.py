#===================================================================================
# All constants and functions utilized in post-processing of DWD populations to
# simulate Milky-Way-like galaxies. Specifically, LISA_FIRE_Lband() creates
# MW-like galaxies of DWDs that are orbiting in the LISA frequency band at present.
#
# Authors: Sarah Thiele & Katelyn Breivik
# Last updated: Oct 2nd, 2021
# Updates: took out snr information
# took out magnitude
# took out ZTF/Tyson
# took out t_evol1 & t_evol2, teff_1, teff_2, etc.
# because we will use the new reduced dat files, the conv files are all correct
# and you don't need to get conv from bpp's for CO+CO and ONe+XX

#===================================================================================
# Imports and Constants:
#===================================================================================

import pandas as pd
import numpy as np
pd.options.mode.chained_assignment = None
from scipy.special import jv
from astropy import constants as const
from astropy import units as u
import astropy.coordinates as coords
from astropy.time import Time
import legwork.utils as utils
import legwork.strain as strain
import legwork.snr as snr
import legwork.psd as lisa
import legwork.evol as evol
from legwork import source
import pdb

G = const.G.value
c = const.c.value  # speed of light in m s^-1
M_sol = const.M_sun.value  # sun's mass in kg
R_sol = const.R_sun.value  # sun's radius in metres
sec_Myr = u.Myr.to('s')  # seconds in a million years
m_kpc = u.kpc.to('m')  # metres in a kiloparsec
L_sol = const.L_sun.value  # solar luminosity in Watts
Z_sun = 0.02  # solar metallicity
sun = coords.get_sun(Time("2021-04-23T00:00:00", scale='utc'))
sun_g = sun.transform_to(coords.Galactocentric)
sun_yGx = sun_g.galcen_distance.to('kpc').value
sun_zGx = sun_g.z.to('kpc').value
M_astro = 7070  # FIRE star particle mass in solar masses
mag_lim = 23  # chosen bolometric magnitude limit

#FIRE = pd.read_hdf('FIRE.h5').sort_values('met')

# Specific to Thiele et al. (2021), here are the used metallicity
# array, the associated binary fractions for each Z value, and the ratios 
# of mass in singles to mass in binaries of the Lband with each specific 
# binary fraction as found using COSMIC's independent samplers
# (See Binary_Fraction_Modeling.ipynb for Tutorials). All values were
# rounded to 4 significant digits except metallicity which used 8:

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

#===================================================================================
# Initial Metallicity/Binary Fraction Functions:
#===================================================================================


def get_Z_from_FeH(FeH, Z_sun=0.02):
    """Converts from FeH to Z under the assumption that
    all stars have the same abundance as the sun
    Parameters
    ----------
    FeH : `array`
        array of Fe/H values to convert
    Z_sun : `float`
        solar metallicity
    Returns
    -------
    Z : `array`
        array of metallicities
    """
    Z = 10**(FeH + np.log10(Z_sun))
    return Z

def get_FeH_from_Z(Z, Z_sun=0.02):
    """Converts from Z to FeH under the assumption that
    all stars have the same abundance as the sun
    Parameters
    ----------
    Z : `array`
        array of metallicities to convert
    Z_sun : `float`
        solar metallicity
    Returns
    -------
    Z : `array`
        array of FeH values
    """
    FeH = np.log10(Z) - np.log10(Z_sun)
    return FeH

def get_binfrac_of_Z(Z):
    '''
    Calculates the theoretical binary fraction as a function
    of metallicity.
    
    Inputs: metallicity Z values
    
    Outputs: binary fraction values
    '''
    FeH = get_FeH_from_Z(Z)
    FeH_low = FeH[np.where(FeH<=-1.0)]
    FeH_high = FeH[np.where(FeH>-1.0)]
    binfrac_low = -0.0648 * FeH_low + 0.3356
    binfrac_high = -0.1977 * FeH_high + 0.2025
    binfrac = np.append(binfrac_low, binfrac_high)
    return binfrac

def get_ratios(binfracs):
    '''
    Calculates the ratio of mass in singles to mass
    in binaries in cosmic for solar metallicity stars.
    These ratios can be used to scale galaxies for a specific
    binary fraction.
    '''
    from cosmic.sample.initialbinarytable import InitialBinaryTable
    from cosmic.sample.sampler import independent
    final_kstar1 = [10, 11, 12]
    final_kstar2 = [10, 11, 12]
    primary_model = 'kroupa01'
    porb_model = 'log_uniform'
    ecc_model = 'uniform'
    SF_start = 13700.0
    SF_duration = 0.0
    met = 0.02
    Size = 1000000
    def InitBinSample(met, size, binfrac):
        InitialBinaries, mass_singles, mass_binaries, n_singles, n_binaries = InitialBinaryTable.sampler('independent', 
                                                                                                         final_kstar1, 
                                                                                                         final_kstar2, 
                                                                                                         binfrac_model=binfrac, 
                                                                                                         primary_model=primary_model, 
                                                                                                         ecc_model=ecc_model, 
                                                                                                         porb_model=porb_model, 
                                                                                                         SF_start=SF_start, 
                                                                                                         SF_duration=SF_duration, 
                                                                                                         met=met, size=size)        
        ratio = mass_singles / mass_binaries
        return ratio

    ratio_05 = InitBinSample(met, size, 0.5)
    ratio_05 = np.round(ratio_05, 2)

    ratios = []
    for binfrac in binfracs:
        ratio = InitBinSample(met, size, binfrac)
        ratios.append(ratio)
    ratios = np.array(ratios)
    ratios = np.round(ratios, 2)

    return ratio_05, ratios


#===================================================================================
# File and Organizing Functions:
#===================================================================================

def getfiles(kstar1, kstar2):
    met_list = [0.0001, 0.00015029, 0.00022588, 0.00033948, 0.00051021, 
                0.00076681, 0.00115245, 0.00173205, 0.00260314, 0.00391233, 
                0.00587993, 0.0088371, 0.01328149, 0.01996108, 0.03]
    
    filename_list = []
    for met in met_list:
        fname = 'dat_kstar1_{}_kstar2_{}_SFstart_13700.0_SFduration_0.0_metallicity_{}.h5'.format(kstar1, kstar2, met)
        filename_list.append(fname)
        
    if kstar1 == '12':
        label='12'
    else:
        label='{}_{}'.format(kstar1, kstar2)
    
    return filename_list, label


def weighted_avg_and_std(values, weights):
    """
    Return the weighted average and standard deviation.

    values, weights -- Numpy ndarrays with the same shape.
    """
    average = np.average(values, weights=weights)
    variance = np.average((values-average)**2, weights=weights)
    return average, np.sqrt(variance)

#===================================================================================
# Lband and Evolution Functions:
#===================================================================================


def beta_(pop):
    '''
    Beta constant from page 8 of Peters(1964) used in the evolution 
    of DWDs due to gravitational waves.
    
    Parameters
    ----------
    pop : `pandas dataframe`
        DF of population which includes component masses in solar
        masses
    Returns
    -------
    beta : `array`
        array of beta values
    '''
    m1 = pop.mass_1 * M_sol
    m2 = pop.mass_2 * M_sol
    beta = 64 / 5 * G ** 3 * m1 * m2 * (m1 + m2) / c ** 5
    return beta

def a_of_t(pop, t):
    '''
    Uses Peters(1964) equation (5.9) for circular binaries to find separation.
    as a function of time.

    Input: the population dataframe from COSMIC. Time t must be in Myr.

    Returns: array of separation at time t in solar radii.
    '''
    t = t * sec_Myr
    beta = beta_(pop)
    a_i = pop.sep * R_sol
    a = (a_i ** 4 - 4 * beta * t) ** (1/4)
    return a / R_sol

def porb_of_a(pop, a):
    '''
    Converts semi-major axis "a" to orbital period using Kepler's equations.

    Input the population dataframe from COSMIC. "a" must be in solar radii and
    an array of the same length as the dateframe pop.

    Returns orbital period in days.
    '''
    a = a * R_sol
    m1 = pop.mass_1 * M_sol
    m2 = pop.mass_2 * M_sol
    P_sqrd = 4 * np.pi ** 2 * a ** 3 / G / (m1 + m2)
    P = np.sqrt(P_sqrd)
    P = P / 3600 / 24  # converts from seconds to days
    return P

def t_of_a(pop, a):
    '''
    Finds time from SRF at which a binary would have a given separation after
    evolving due to gw radiation. (Re-arrangement of a_of_t(pop, t)).

    "a" must be in solar radii.

    Returns time in Myr.
    '''
    beta = beta_(pop)
    a_i = pop.sep * R_sol
    a = a * R_sol
    t = (a_i ** 4 - a ** 4) / 4 / beta
    t = t / sec_Myr
    return t

def t_merge(pop):
    '''
    Uses Peters(1964) equation (5.10) to determine the merger time of a circular
    DWD binary from time of SRF.

    Returns time in Myr.
    '''
    a_0 = pop.sep * R_sol
    beta = beta_(pop)
    T = a_0 ** 4 / 4 / beta
    T / sec_Myr
    return T

def a_of_RLOF(set):
    '''
    Finds separation when secondary overflows its
    Roche Lobe. Returns "a" in solar radii.
    
    Taken from Eq. 23 in "Binary evolution in a nutshell" 
    by Marc van der Sluys, which is an approximation of a fit
    done of Roche-lobe radius by Eggleton (1983).
    '''
    m1 = set.mass_1
    m2 = set.mass_2
    R2 = set.rad_2
    q = m2 / m1
    num = 0.49 * q ** (2/3)
    denom = 0.6 * q ** (2/3) + np.log(1 + q ** (1/3))
    a = denom * R2 / num
    return a

def random_sphere(R, num):
    '''
    Generates "num" number of random points within a
    sphere of radius R. It picks random x, y, z values
    within a cube and discards it if it's outside the
    sphere.

    Inputs: Radius in kpc, num is an integer

    Outputs: X, Y, Z arrays of length num
    '''
    X = []
    Y = []
    Z = []
    while len(X) < num:
        x = np.random.uniform(-R, R)
        y = np.random.uniform(-R, R)
        z = np.random.uniform(-R, R)
        r = np.sqrt(x ** 2 + y ** 2 + z ** 2)
        if r > R:
            continue
        if r <= R:
            X.append(x)
            Y.append(y)
            Z.append(z)
    X = np.array(X)
    Y = np.array(Y)
    Z = np.array(Z)
    return X, Y, Z

def rad_WD(M):
    '''
    Calculates the radius of a WD as a function of mass M in solar masses.
    Taken from Eq. 91 in Hurley et al. (2000), from Eq. 17 in Tout et al. (1997)
    
    Input an array M of mass in solar masses
    
    Outputs the radius of the WD in solar radii
    '''
    M_ch = 1.44
    R_NS = 1.4e-5*np.ones(len(M))
    A = 0.0115 * np.sqrt((M_ch/M)**(2/3) - (M/M_ch)**(2/3))
    rad = np.max(np.array([R_NS, A]), axis=0)
    return rad

def A_val(kstar): 
    '''
    Returns the baryon number for each WD type
    to be used in (modified) Mestel cooling in
    WD_Cooling function. 4 corresponds to He,
    15 to CO and 17 to ONe.
    '''
    A_list = np.zeros(len(kstar)) 
    A_list[np.where(kstar == 10)] = 4 
    A_list[np.where(kstar == 11)] = 15 
    A_list[np.where(kstar == 12)] = 17 
    return A_list 

def WD_Cooling(data, i):
    '''
    Modified Mestel Cooling as specified in Hurley (2003), 
    equation 1. Returns the evolved present-day luminosity
    of each binary component.
    
    Returns luminosity in Watts.
    '''
    def L_WD(A, b, x, M, Z, t):
        return (b * M * Z ** (0.4)) / (A * (t + 0.1)) ** x
    
    def b_old(A):
        b = 300 * (9000 * A) ** (5.3)
        return b
    
    x_old = 6.48
    b_young = 300
    x_young = 1.18

    A1 = A_val(data.kstar_1.values)
    A2 = A_val(data.kstar_2.values)
    M1 = data.mass_1.values
    M2 = data.mass_2.values
    Z = met_arr[i+1] * np.ones(len(data))    
    t1 = data.t_evol_1.values
    t2 = data.t_evol_2.values
    
    L1_young = L_WD(A1, b_young, x_young, M1, Z, t1)
    L1 = L_WD(A1, b_old(A1), x_old, M1, Z, t1)
    L2_young = L_WD(A2, b_young, x_young, M2, Z, t2)
    L2 = L_WD(A2, b_old(A2), x_old, M2, Z, t2)
    
    L1[t1 < 9000.0] = L1_young[t1 < 9000.0]
    L2[t2 < 9000.0] = L2_young[t2 < 9000.0]
    L1 = L1 * L_sol
    L2 = L2 * L_sol
    return L1, L2

def T_eff(L1, L2, data):
    '''
    Calculates effective temperature using Stefan-Boltzmann
    Law.
    '''
    sigma = 5.67037e-8 # Stefan-Boltzmann Constant
    r1 = data.rad_1.values * R_sol
    r2 = data.rad_2.values * R_sol
    T1 = (L1 / (4 * np.pi * sigma * r1 ** 2)) ** (1 / 4)
    T2 = (L2 / (4 * np.pi * sigma * r2 ** 2)) ** (1 / 4)
    return T1, T2

def mag_bol(data, i):
    L1, L2 = WD_Cooling(data, i)
    d = data.dist_sun.values * 1000
    M1_bol = 4.8 - 2.5 * np.log10(L1 / L_sol)
    m1_bol = M1_bol + 5 * np.log10(d / 10)
    M2_bol = 4.8 - 2.5 * np.log10(L2 / L_sol)
    m2_bol = M2_bol + 5 * np.log10(d / 10)
    m_tot = -2.5 * np.log10(10 ** (-0.4 * m1_bol) + 10 ** (-0.4 * m2_bol))
    return m_tot

def LISA_calcs(LISA_band):
    '''
    get h_0, chirp mass and chirp values for LISA-band systems.
    '''
    f_orb = LISA_band.f_gw.values / 2 * u.s**(-1)
    ecc = np.zeros(len(LISA_band))
    m1 = LISA_band.mass_1.values * u.M_sun
    m2 = LISA_band.mass_2.values * u.M_sun
    mc = utils.chirp_mass(m1, m2)
    t_obs = 4*u.yr
    dist = LISA_band.dist_sun.values * u.kpc   
    sources = source.Source(m_1=m1, m_2=m2, ecc=ecc, dist=dist, f_orb=f_orb,
                            gw_lum_tol=0.05, stat_tol=1e-2, interpolate_g=True)
    h_0 = sources.get_h_0_n(harmonics=[2]).reshape(len(LISA_band))
    chirps = utils.fn_dot(mc, f_orb, ecc, 2)
    
    LISA_band['h_0'] = h_0
    LISA_band['fdot'] = chirps.value
    LISA_band['']
    
    return LISA_band

def evolve(pop_init):
    '''
    Evolve an initial population of binary WD's using
    GW radiation.
    '''
    t_evol = pop_init.age * 1000 - pop_init.tphys
    sep_f = a_of_t(pop_init, t_evol)
    porb_f = porb_of_a(pop_init, sep_f)
    f_gw = 2 / (porb_f * 24 * 3600)
    pop_init['t_evol'] = t_evol
    pop_init['sep_f'] = sep_f
    pop_init['porb_f'] = porb_f
    pop_init['f_gw'] = f_gw
    return pop_init

def position(pop_init):
    '''
    Assigning random microchanges to positions to
    give each system a unique position for identical
    FIRE star particles
    '''
    R_list = pop_init.kern_len.values
    xGx = pop_init.xGx.values.copy()
    yGx = pop_init.yGx.values.copy()
    zGx = pop_init.zGx.values.copy()
    x, y, z = random_sphere(1.0, len(R_list))
    X = xGx + (x * R_list)
    Y = yGx + (y * R_list)
    Z = zGx + (z * R_list)
    pop_init['X'] = X
    pop_init['Y'] = Y
    pop_init['Z'] = Z
    pop_init['dist_sun'] = (X ** 2 + (Y - sun_yGx) ** 2 + (Z - sun_zGx) ** 2) ** (1/2)   
    return pop_init
  
def merging_pop(pop_init):
    t_m = t_merge(pop_init)
    pop_init['t_delay'] = t_m + pop_init.tphys
    pop_merge = pop_init.loc[pop_init.t_delay <= pop_init.age * 1000]
    pop_init = pop_init.loc[pop_init.t_delay >= pop_init.age * 1000]
    return pop_init, pop_merge

def RLOF_pop(pop_init):
    a_RLOF = a_of_RLOF(pop_init)
    t_RLOF = t_of_a(pop_init, a_RLOF)
    pop_init['t_RLOF'] = t_RLOF
    pop_RLOF = pop_init.loc[t_RLOF + pop_init.tphys <= pop_init.age * 1000]
    pop_init = pop_init.loc[t_RLOF + pop_init.tphys >= pop_init.age * 1000]
    return pop_init, pop_RLOF
    
#===================================================================================
# Functions to create LISA galaxies:
#===================================================================================


def create_population(pop_init, i, label, ratio, binfrac, interfile):
    pop_init[['bin_num', 'FIRE_index']] = pop_init[['bin_num', 'FIRE_index']].astype('int64')
    if interfile == True:
        pop_init[['bin_num', 'FIRE_index']].to_hdf('Lband_{}_{}_{}_inter.hdf'.format(label,
                                                                                     met_arr[i+1],
                                                                                     binfrac),
                                                   key='pop_init', format='t', append=True)    
    # Now that we've obtained an initial population, we make data cuts
    # of systems who wouldn't form in time for their FIRE age, or would
    # merge or overflow their Roche Lobe before present day.
    pop_init = pop_init.loc[pop_init.tphys <= pop_init.age * 1000]
    if interfile == True:
        pop_init[['bin_num', 'FIRE_index']].to_hdf('Lband_{}_{}_{}_inter.hdf'.format(label, 
                                                                                     met_arr[i+1], 
                                                                                     binfrac), 
                                                    key='pop_age', format='t', append=True)
    
    pop_init, pop_merge = merging_pop(pop_init)
    if interfile == True:
        pop_merge[['bin_num', 'FIRE_index']].to_hdf('Lband_{}_{}_{}_inter.hdf'.format(label,
                                                                                      met_arr[i+1], 
                                                                                      binfrac), 
                                                    key='pop_merge', format='t', append=True)    
    
        pop_init[['bin_num', 'FIRE_index']].to_hdf('Lband_{}_{}_{}_inter.hdf'.format(label, 
                                                                                     met_arr[i+1], 
                                                                                     binfrac), 
                                key='pop_nm', format='t', append=True)
    
    pop_merge = pd.DataFrame()
    pop_init, pop_RLOF = RLOF_pop(pop_init)
    
    if interfile == True:
        pop_RLOF[['bin_num','FIRE_index']].to_hdf('Lband_{}_{}_{}_inter.hdf'.format(label,
                                                                                           met_arr[i+1], 
                                                                                           binfrac), 
                                                  key='pop_RLOF', format='t', append=True)
    
        pop_init[['bin_num', 'FIRE_index']].to_hdf('Lband_{}_{}_{}_inter.hdf'.format(label, 
                                                                                            met_arr[i+1], 
                                                                                            binfrac), 
                                key='pop_nRLOF', format='t', append=True)
    pop_RLOF = pd.DataFrame()
    
    # We now have a final population which we can evolve
    # using GW radiation
    pop_init = evolve(pop_init)
    
    # Assigning random microchanges to positions to
    # give each system a unique position for identical
    # FIRE star particles
    pop_init = position(pop_init)
    
    if interfile == True:
        pop_init[['bin_num', 'FIRE_index', 'X', 'Y', 'Z']].to_hdf('Lband_{}_{}_{}_inter.hdf'.format(label, 
                                                                                                    met_arr[i+1], 
                                                                                                    binfrac), 
                                                                  key='pop_f', format='t', append=True)    
    if binfrac == 0.5:
        binfrac_write = 0.5
    else:
        binfrac_write = 'variable'
    
    # Assigning weights to population to be used for histograms.
    # This creates an extra columns which states how many times
    # a given system was sampled from the cosmic-pop conv df.
    pop_init = pop_init.join(pop_init.groupby('bin_num')['bin_num'].size(), 
                             on='bin_num', rsuffix='_pw')
    
    # Systems detectable by LISA will be in the frequency band
    # between f_gw's 0.01mHz and 1Hz.
    LISA_band = pop_init.loc[(pop_init.f_gw >= 1e-4)]
    if len(LISA_band) == 0:
        print('No LISA sources')
        return
    else:
        pop_init = pd.DataFrame()
        LISA_band = LISA_band.join(LISA_band.groupby('bin_num')['bin_num'].size(), 
                                   on='bin_num', rsuffix='_Lw')
        print('got LISA band and added weight column')
    
        # LISA GW calculations:
        #LISA_band = LISA_calcs(LISA_band)
        #print('finished LISA band calculations')

        # Output to hdf files
        LISA_band.to_hdf('Lband_{}_{}_{}.hdf'.format(label, 
                                                     met_arr[i+1], binfrac), 
                         key='Lband', format='t', append=True)
    
        return
    
    
    
def make_galaxy(path, filename, i, label, ratio, binfrac, interfile, fire_path):
    FIRE = pd.read_hdf(fire_path+'FIRE.h5').sort_values('met')

    rand_seed = np.random.randint(0, 100, 1)
    np.random.seed(rand_seed)
    
    rand_seed = pd.DataFrame(rand_seed)
    rand_seed.to_hdf('Lband_{}_{}_{}.hdf'.format(label, met_arr[i+1], binfrac), 
                     key='rand_seed')

    # Choose metallicity bin
    met_start = met_arr[i] / Z_sun
    met_end = met_arr[i+1] / Z_sun
    
    # Calculating the formation time of each component:
    conv = pd.read_hdf(path+filename, key='conv')
    
    # Re-writing the radii of each component since the conv df 
    # doesn't log the WD radius properly
    conv['rad_1'] = rad_WD(conv.mass_1.values)
    conv['rad_2'] = rad_WD(conv.mass_2.values)    
    
    # Use ratio to scale to astrophysical pop w/ specific binary frac.
    try:
        mass_binaries = pd.read_hdf(path+filename, key='mass_stars').iloc[-1]
    except:
        print('m_binaries key')
        mass_binaries = pd.read_hdf(path+filename, key='mass_binaries').iloc[-1]
    mass_total = (1 + ratio) * mass_binaries
    
    mass_total.to_hdf('Lband_{}_{}_{}.hdf'.format(label, met_arr[i+1], 
                                                  binfrac), key='mass_total')
    DWD_per_mass = len(conv) / mass_total
    N_astro = DWD_per_mass * M_astro  # num of binaries per star particle
    
    # Choose FIRE bin based on metallicity
    FIRE['FIRE_index'] = FIRE.index
    if met_end * Z_sun == met_arr[-1]:
        FIRE_bin = FIRE.loc[FIRE.met >= met_start]
    else:
        FIRE_bin = FIRE.loc[(FIRE.met >= met_start)&(FIRE.met <= met_end)]
    FIRE = []
    
    # We sample by the integer number of systems per star particle,
    # as well as a probabilistic approach for the fractional component
    # of N_astro:
    N_astro_dec = N_astro % 1
    p_DWD = np.random.rand(len(FIRE_bin))
    N_sample_dec = np.zeros(len(FIRE_bin))
    N_sample_dec[p_DWD <= N_astro_dec.values] = 1.0
    num_sample_dec = int(N_sample_dec.sum())
    print('we will sample {} stars from the decimal portion'.format(num_sample_dec))
    sample_dec = pd.DataFrame.sample(conv, num_sample_dec, replace=True)
    FIRE_bin2 = FIRE_bin.loc[N_sample_dec == 1.0]
    
    params_list = ['bin_num', 'mass_1', 'mass_2', 'kstar_1', 'kstar_2', 'porb', 'sep', 
            'met', 'age', 'tphys', 'rad_1', 'rad_2', 'kern_len', 'xGx', 'yGx', 'zGx', 
                   'FIRE_index']#, 'CEsep', 'CEtime', 'RLOFsep', 'RLOFtime']
    
    pop_init = pd.concat([sample_dec.reset_index(), FIRE_bin2.reset_index()], axis=1)
    sample_dec = pd.DataFrame()
    FIRE_bin2 = pd.DataFrame()
    create_population(pop_init[params_list], i, label, ratio, binfrac, interfile)

    N_sample_int = int(N_astro) * len(FIRE_bin)
    print('we will sample {} stars from the integer portion'.format(N_sample_int))

    print('getting FIRE bin')
    FIRE_repeat = pd.DataFrame(np.repeat(FIRE_bin.values, int(N_astro), axis=0))
    FIRE_repeat.columns = FIRE_bin.columns
    FIRE_bin = pd.DataFrame()
    
    if N_sample_int < 1e6:
        sample_int = pd.DataFrame.sample(conv, N_sample_int, replace=True)
        pop_init_int = pd.concat([sample_int.reset_index(), 
                              FIRE_repeat.reset_index()], axis=1)
        N = len(sample_int)
        sample_int = pd.DataFrame()
        FIRE_repeat = pd.DataFrame()
        create_population(pop_init_int[params_list], i, label, ratio, binfrac, interfile)

    elif N_sample_int > 1e6:
        print('looping the integer population')
        N = 0
        j = 0
        jlast = int(1e6)
        while j < N_sample_int:
            print('j: ', j)
            print('jlast: ', jlast)
            print('sampling {} systems'.format(int(jlast - j)))
            sample_int = pd.DataFrame.sample(conv, int(jlast-j), replace=True)
            N += len(sample_int)
            pop_init_int = pd.concat([sample_int.reset_index(), 
                                      FIRE_repeat.iloc[j:jlast].reset_index()], axis=1)
            create_population(pop_init_int[params_list], i, label, ratio, binfrac, interfile)
            j += 1e6
            j = int(j)
            jlast += 1e6
            jlast = int(jlast)
            if jlast > N_sample_int:
                jlast = N_sample_int
    if N != N_sample_int:
        print('loop is incorrect')

    FIRE_repeat = pd.DataFrame()
    sample_int = pd.DataFrame()
    
    return

def Lband_files(kstar1, kstar2, var=True):
    met_list = [0.0001, 0.00015029, 0.00022588, 0.00033948, 0.00051021, 
                0.00076681, 0.00115245, 0.00173205, 0.00260314, 0.00391233, 
                0.00587993, 0.0088371, 0.01328149, 0.01996108, 0.03]
    if var:
        var_list = [0.4847, 0.4732, 0.4618, 0.4503, 0.4388, 
                    0.4274, 0.4159, 0.4044, 0.3776, 0.3426, 
                    0.3076, 0.2726, 0.2376, 0.2027, 0.1677]
        
        if kstar1 != '12':
            files = ['Lband_{}_{}_{}_{}.hdf'.format(kstar1, kstar2, met, var) for var, met in zip(var_list, met_list)]
        else:
            files = ['Lband_{}_{}_{}.hdf'.format(kstar1, met, var) for var, met in zip(var_list, met_list)]
    else:
        if kstar1 != '12':
            files = ['Lband_{}_{}_{}_0.5.hdf'.format(kstar1, kstar2, met) for met in met_list]
        else:
            files = ['Lband_{}_{}_0.5.hdf'.format(kstar1, met) for met in met_list]
    
    return files

def galaxy_files(kstar1, kstar2, var=True):
    met_list = [0.0001, 0.00015029, 0.00022588, 0.00033948, 0.00051021, 
                0.00076681, 0.00115245, 0.00173205, 0.00260314, 0.00391233, 
                0.00587993, 0.0088371, 0.01328149, 0.01996108, 0.03]
    if var:
        var_list = [0.4847, 0.4732, 0.4618, 0.4503, 0.4388, 
                    0.4274, 0.4159, 0.4044, 0.3776, 0.3426, 
                    0.3076, 0.2726, 0.2376, 0.2027, 0.1677]
        
        if kstar1 != '12':
            files = ['final_galaxy_{}_{}_{}_{}.hdf'.format(kstar1, kstar2, met, var) for var, met in zip(var_list, met_list)]
        else:
            files = ['final_galaxy_{}_{}_{}.hdf'.format(kstar1, met, var) for var, met in zip(var_list, met_list)]
    else:
        if kstar1 != '12':
            files = ['final_galaxy_{}_{}_{}_0.5.hdf'.format(kstar1, kstar2, met) for met in met_list]
        else:
            files = ['final_galaxy_{}_{}_0.5.hdf'.format(kstar1, met) for met in met_list]
    
    return files
