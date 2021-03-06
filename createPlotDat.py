#=========================================================================
# This script creates LISA DWD galaxies across 15 metallicity bins,
# incorporating the metallicity-dependent binary fraction as
# discussed in Thiele et al. (2021). It is the ongoing,
# unstable version. Older versions can be found at (insert link).
#
# Author: Sarah Thiele
# Last updated: Oct 2nd, 2021
# Last changes to create.py: take out ZTF/Tyson files, add in option to not write-out to 
# interfiles. New functions module (funcs_v1)
#=========================================================================

import numpy as np
from astropy import constants as const
from astropy import units as u
import astropy.coordinates as coords
from astropy.time import Time
import argparse
import postproc as pp
import utils


# Set constants:
# LEGWORK uses astropy units so we do also for consistency
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
parser.add_argument('--path', default='./', help='path to COSMIC dat files')
parser.add_argument('--lband-path', default='./', help='path to save LISA band DWD data')
parser.add_argument('--plotdat-path', default='./', help='path to save plotting data')
args = parser.parse_args()

pp.get_formeff(args.path, args.lband_path, args.plotdat_path, getfrom='dat')

pp.get_interactionsep(args.path, args.lband_path, args.plotdat_path, verbose=False)
pp.get_numLISA(args.lband_path, args.plotdat_path, Lbandfile='new', FIREmin=0.00015, FIREmax=13.346, Z_sun=0.02)

pp.get_resolvedDWDs(args.lband_path, args.plotdat_path, var=True, window=1000)
pp.get_resolvedDWDs(args.lband_path, args.plotdat_path, var=False, window=1000)
