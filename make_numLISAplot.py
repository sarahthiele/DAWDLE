from funcs_v1
import matplotlib.pyplot as plt
import legwork.visualisation as vis
from matplotlib.colors import TwoSlopeNorm
import seaborn as sb
from matplotlib import rcParams

rcParams['font.family'] = 'serif'
rcParams['font.size'] = 14
rcParams['mathtext.default'] = 'regular'

def make_numLISAplot(numsFZ, numsF50):
    num = 30
    met_bins = np.logspace(np.log10(FIRE.met.min()), np.log10(FIRE.met.max()), num)*Z_sun
    met_bins

    nums = pd.read_hdf('numLISA_30bins.hdf', key='data')
    nums05 = pd.read_hdf('numLISA_30bins_05.hdf', key='data')

    Henums = numsFZ.He.values
    COHenums = numsFZ.COHe.values
    COnums = numsFZ.CO.values
    ONenums = numsFZ.ONe.values

    Henums05 = numsF50.He.values
    COHenums05 = numsF50.COHe.values
    COnums05 = numsF50.CO.values
    ONenums05 = numsF50.ONe.values

    fig, ax = plt.subplots(1, 4, figsize=(24, 6))

    ax[0].plot(np.log10(met_bins[1:]/Z_sun), Henums/1e5, drawstyle='steps-mid', 
               color='xkcd:tomato red', lw=3, label='f$_b$(Z)')
    ax[0].plot(np.log10(met_bins[1:]/Z_sun), Henums05/1e5, 
               drawstyle='steps-mid', color='xkcd:tomato red', ls='--', lw=3, label='f$_b$=0.5')
    ax[0].text(0.05, 0.9, 'He + He', fontsize=20, transform=ax[0].transAxes)

    ax[1].plot(np.log10(met_bins[1:]/Z_sun), COHenums/1e5, drawstyle='steps-mid', 
               color='xkcd:blurple', lw=3, label='f$_b$(Z)')
    ax[1].plot(np.log10(met_bins[1:]/Z_sun), COHenums05/1e5, drawstyle='steps-mid', 
               color='xkcd:blurple', ls='--', lw=3, label='f$_b$=0.5')
    ax[1].text(0.05, 0.9, 'CO + He', fontsize=20, transform=ax[1].transAxes)

    ax[2].plot(np.log10(met_bins[1:]/Z_sun), COnums/1e5, drawstyle='steps-mid', 
               color='xkcd:pink', lw=3, label='f$_b$(Z)')
    ax[2].plot(np.log10(met_bins[1:]/Z_sun), COnums05/1e5, drawstyle='steps-mid', 
               color='xkcd:pink', ls='--', lw=3, label='f$_b$=0.5')
    ax[2].text(0.05, 0.9, 'CO + CO', fontsize=20, transform=ax[2].transAxes)

    ax[3].plot(np.log10(met_bins[1:]/Z_sun), ONenums/1e5, drawstyle='steps-mid', 
               color='xkcd:light blue', lw=3, label='f$_b$(Z)')
    ax[3].plot(np.log10(met_bins[1:]/Z_sun), ONenums05/1e5, drawstyle='steps-mid',
               color='xkcd:light blue', ls='--', lw=3, label='f$_b$=0.5')
    ax[3].text(0.05, 0.9, 'ONe + X', fontsize=20, transform=ax[3].transAxes)

    for i in range(4):
        #ax[i].set_yscale('log')
        #ax[i].set_ylim(10, 2.5e6)
        #ax[i].grid(which='both', zorder=0, alpha=0.2)
        ax[i].set_xlabel('Log$_{10}$(Z/Z$_\odot$)')
        ax[i].set_xticks([-3, -2, -1, 0, 1.])
        ax[i].legend(loc='lower left', bbox_to_anchor= (-0.02, 1.01), ncol=2, 
                     borderaxespad=0, frameon=False, fontsize=21)

    plt.subplots_adjust(wspace=0.2)   
    ax[0].set_ylabel(r'N$_{\rm{DWD}}$(f$_{\rm{GW}} \geq 10^{-4} \rm{HZ}$) (10$^5$)')
    ax[0].set_yticks(np.arange(0, 2.5, 0.5));
    ax[2].set_yticks(np.arange(0, 3.5, 0.5));
    plt.show(block=False)
    
    return 
