#=========================================================================
# This script creates LISA DWD galaxies across 13 metallicity bins,
# incorporating the metallicity-dependent binary fraction as
# discussed in Thiele et al. (2021). It is the ongoing,
# unstable version. Older versions can be found at (insert link).
#
# Author: Sarah Thiele
# Last updated: April 20th, 2021
#=========================================================================

import pandas as pd
import numpy as np
from cosmic import MC_samp
pd.options.mode.chained_assignment = None
from scipy.special import jv
from astropy import constants as const
from astropy import units as u
import astropy.coordinates as coords
from astropy.time import Time
import argparse
from functions import *
import legwork.utils as utils
import legwork.strain as strain
import legwork.snr as snr
import legwork.lisa as lisa
import legwork.evol as evol
from legwork import source

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
parser = argparse.ArgumentParser()
parser.add_argument("--DWD", default="He_He", type=str)
args = parser.parse_args()

path = '/mnt/raid-cita/sthiele/DWD_Met_Grid/'
FIRE = pd.read_hdf(path + 'FIRE.h5').sort_values('met')

met_arr = np.logspace(np.log10(1e-4), np.log10(0.03), 15)
met_arr = np.round(met_arr, 8)
met_arr = np.append(0.0, met_arr)

binfracs = np.array([0.4847, 0.4732, 0.4618, 0.4503, 0.4388, 0.4274, 0.4159, 0.4044, 
		     0.3776, 0.3426, 0.3076, 0.2726, 0.2376, 0.2027, 0.1677])

ratios = np.array([0.68, 0.71, 0.74, 0.78, 0.82, 0.86, 0.9, 
		   0.94, 1.05, 1.22, 1.44, 1.7 , 2.05, 2.51, 3.17])

ratio_05 = 0.64

ZTF = True
Tyson = True

# Run Code:
if args.DWD == 'He_He':
    fname, label = getfiles_He_He(path)
elif args.DWD == 'CO_He':
        fname, label = getfiles_CO_He(path)
elif args.DWD == 'CO_CO':
    fname, label = getfiles_CO_CO(path)
elif args.DWD == 'ONe':
    fname, label = getfiles_ONe(path)

i = 0
for f in fname:
    ratio = ratios[i] 
    binfrac = binfracs[i]
    LISA_FIRE_galaxy(f, i, label, ratio, binfrac, ZTF, Tyson)
    LISA_FIRE_galaxy(f, i, label, ratio_05, 0.5, ZTF, Tyson)
    i += 1


