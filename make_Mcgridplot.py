from funcs_v1 import *
obs_sec = 4 * u.yr.to('s')
obs_hz = 1 / obs_sec

def make_Mc_fgw_plot(pathtoLband, model):
    if model == 'FZold':
        He = pd.DataFrame()
        for f in galaxy_files_10_10_var():
            He = He.append(pd.read_hdf(pathtoLband + f, key='Lband'))

        CO = pd.DataFrame()
        for f in galaxy_files_11_11_var():
            CO = CO.append(pd.read_hdf(pathtoLband + f, key='Lband'))

        COHe = pd.DataFrame()
        for f in galaxy_files_11_10_var():
            COHe = COHe.append(pd.read_hdf(pathtoLband + f, key='Lband'))

        ONe = pd.DataFrame()
        for f in galaxy_files_12_var():
            ONe = ONe.append(pd.read_hdf(pathtoLband + f, key='Lband'))
            
    elif model == 'FZnew':
        He = pd.DataFrame()
        for f in Lband_files_10_10_var():
            He = He.append(pd.read_hdf(pathtoLband + f, key='Lband'))

        CO = pd.DataFrame()
        for f in Lband_files_11_11_var():
            CO = CO.append(pd.read_hdf(pathtoLband + f, key='Lband'))

        COHe = pd.DataFrame()
        for f in Lband_files_11_10_var():
            COHe = COHe.append(pd.read_hdf(pathtoLband + f, key='Lband'))

        ONe = pd.DataFrame()
        for f in Lband_files_12_var():
            ONe = ONe.append(pd.read_hdf(pathtoLband + f, key='Lband'))
            
    if model == 'F50old':
        He = pd.DataFrame()
        for f in galaxy_files_10_10_05():
            He = He.append(pd.read_hdf(pathtoLband + f, key='Lband'))

        CO = pd.DataFrame()
        for f in galaxy_files_11_11_05():
            CO = CO.append(pd.read_hdf(pathtoLband + f, key='Lband'))

        COHe = pd.DataFrame()
        for f in galaxy_files_11_10_05():
            COHe = COHe.append(pd.read_hdf(pathtoLband + f, key='Lband'))

        ONe = pd.DataFrame()
        for f in galaxy_files_12_05():
            ONe = ONe.append(pd.read_hdf(pathtoLband + f, key='Lband'))
            
    elif model == 'F50new':
        He = pd.DataFrame()
        for f in Lband_files_10_10_05():
            He = He.append(pd.read_hdf(pathtoLband + f, key='Lband'))

        CO = pd.DataFrame()
        for f in Lband_files_11_11_05():
            CO = CO.append(pd.read_hdf(pathtoLband + f, key='Lband'))

        COHe = pd.DataFrame()
        for f in Lband_files_11_10_05():
            COHe = COHe.append(pd.read_hdf(pathtoLband + f, key='Lband'))

        ONe = pd.DataFrame()
        for f in Lband_files_12_05():
            ONe = ONe.append(pd.read_hdf(pathtoLband + f, key='Lband'))

    Heplot = He.loc[(He.fdot>=obs_hz)&(He.snr>7)] #[::100]
    COHeplot = COHe.loc[(COHe.fdot>=obs_hz)&(COHe.snr>7)] #[::1000]
    COplot = CO.loc[(CO.fdot>=obs_hz)&(CO.snr>7)] #[::100]
    ONeplot = ONe.loc[(ONe.fdot>=obs_hz)&(ONe.snr>7)] #[::100]

    fig, ax = plt.subplots(4, 3, figsize=(20,16))
    levels = [0.1, 0.3, 0.6, 0.9]
    colors = ['#80afd6', '#2b5d87', '#4288c2', '#17334a']

    ax[0,0].scatter(y=utils.chirp_mass(Heplot.mass_1.values*u.M_sun, 
                                  Heplot.mass_2.values*u.M_sun),
               x=np.log10(Heplot.f_gw.values), color='xkcd:light grey', zorder=0.)

    sb.kdeplot(y=utils.chirp_mass(He.loc[He.met*Z_sun<=met_arr[1]].mass_1.values*u.M_sun, 
                                  He.loc[He.met*Z_sun<=met_arr[1]].mass_2.values*u.M_sun)[::10],
               x=np.log10(He.loc[He.met*Z_sun<=met_arr[1]].f_gw.values)[::10], levels=levels,fill=False, 
               ax=ax[0,0], color=colors[0], zorder=3, linewidths=2.5)

    ax[0,1].scatter(y=utils.chirp_mass(Heplot.mass_1.values*u.M_sun, 
                                  Heplot.mass_2.values*u.M_sun),
               x=np.log10(Heplot.f_gw.values), color='xkcd:light grey', zorder=0.)

    sb.kdeplot(y=utils.chirp_mass(He.loc[(He.met*Z_sun>=met_arr[7])&(He.met*Z_sun<=met_arr[8])].mass_1.values*u.M_sun, 
                                  He.loc[(He.met*Z_sun>=met_arr[7])&(He.met*Z_sun<=met_arr[8])].mass_2.values*u.M_sun)[::10],
               x=np.log10(He.loc[(He.met*Z_sun>=met_arr[7])&(He.met*Z_sun<=met_arr[8])].f_gw.values)[::10], levels=levels,fill=False, 
               ax=ax[0,1], color=colors[1], zorder=3, linewidths=2.5)

    ax[0,2].scatter(y=utils.chirp_mass(Heplot.mass_1.values*u.M_sun, 
                                  Heplot.mass_2.values*u.M_sun),
               x=np.log10(Heplot.f_gw.values), color='xkcd:light grey', zorder=0.)

    sb.kdeplot(y=utils.chirp_mass(He.loc[(He.met*Z_sun>=met_arr[-2])].mass_1.values*u.M_sun, 
                                  He.loc[(He.met*Z_sun>=met_arr[-2])].mass_2.values*u.M_sun)[::100],
               x=np.log10(He.loc[(He.met*Z_sun>=met_arr[-2])].f_gw.values)[::100], levels=levels,fill=False, 
               ax=ax[0,2], color=colors[3], zorder=3, linewidths=2.5)

    ax[1,0].scatter(y=utils.chirp_mass(COHeplot.mass_1.values*u.M_sun, 
                                  COHeplot.mass_2.values*u.M_sun),
               x=np.log10(COHeplot.f_gw.values), color='xkcd:light grey', zorder=0.)

    sb.kdeplot(y=utils.chirp_mass(COHe.loc[COHe.met*Z_sun<=met_arr[1]].mass_1.values*u.M_sun, 
                                  COHe.loc[COHe.met*Z_sun<=met_arr[1]].mass_2.values*u.M_sun)[::10],
               x=np.log10(COHe.loc[COHe.met*Z_sun<=met_arr[1]].f_gw.values)[::10], levels=levels,fill=False, 
               ax=ax[1,0], color=colors[0], zorder=3, linewidths=2.5)

    ax[1,1].scatter(y=utils.chirp_mass(COHeplot.mass_1.values*u.M_sun, 
                                  COHeplot.mass_2.values*u.M_sun),
               x=np.log10(COHeplot.f_gw.values), color='xkcd:light grey', zorder=0.)

    sb.kdeplot(y=utils.chirp_mass(COHe.loc[(COHe.met*Z_sun>=met_arr[7])&(COHe.met*Z_sun<=met_arr[8])].mass_1.values*u.M_sun, 
                                  COHe.loc[(COHe.met*Z_sun>=met_arr[7])&(COHe.met*Z_sun<=met_arr[8])].mass_2.values*u.M_sun)[::100],
               x=np.log10(COHe.loc[(COHe.met*Z_sun>=met_arr[7])&(COHe.met*Z_sun<=met_arr[8])].f_gw.values)[::100], levels=levels,fill=False, 
               ax=ax[1,1], color=colors[1], zorder=3, linewidths=2.5)

    ax[1,2].scatter(y=utils.chirp_mass(COHeplot.mass_1.values*u.M_sun, 
                                  COHeplot.mass_2.values*u.M_sun),
               x=np.log10(COHeplot.f_gw.values), color='xkcd:light grey', zorder=0.)

    sb.kdeplot(y=utils.chirp_mass(COHe.loc[(COHe.met*Z_sun>=met_arr[-2])].mass_1.values*u.M_sun, 
                                  COHe.loc[(COHe.met*Z_sun>=met_arr[-2])].mass_2.values*u.M_sun)[::1000],
               x=np.log10(COHe.loc[(COHe.met*Z_sun>=met_arr[-2])].f_gw.values)[::1000], levels=levels,fill=False, 
               ax=ax[1,2], color=colors[3], zorder=3, linewidths=2.5)

    ax[2,0].scatter(y=utils.chirp_mass(COplot.mass_1.values*u.M_sun, 
                                  COplot.mass_2.values*u.M_sun),
               x=np.log10(COplot.f_gw.values), color='xkcd:light grey', zorder=0.)

    sb.kdeplot(y=utils.chirp_mass(CO.loc[CO.met*Z_sun<=met_arr[1]].mass_1.values*u.M_sun, 
                                  CO.loc[CO.met*Z_sun<=met_arr[1]].mass_2.values*u.M_sun)[::1],
               x=np.log10(CO.loc[CO.met*Z_sun<=met_arr[1]].f_gw.values)[::1], levels=levels,fill=False, 
               ax=ax[2,0], color=colors[0], zorder=3, linewidths=2.5)

    ax[2,1].scatter(y=utils.chirp_mass(COplot.mass_1.values*u.M_sun, 
                                  COplot.mass_2.values*u.M_sun),
               x=np.log10(COplot.f_gw.values), color='xkcd:light grey', zorder=0.)

    sb.kdeplot(y=utils.chirp_mass(CO.loc[(CO.met*Z_sun>=met_arr[7])&(CO.met*Z_sun<=met_arr[8])].mass_1.values*u.M_sun, 
                                  CO.loc[(CO.met*Z_sun>=met_arr[7])&(CO.met*Z_sun<=met_arr[8])].mass_2.values*u.M_sun)[::1],
               x=np.log10(CO.loc[(CO.met*Z_sun>=met_arr[7])&(CO.met*Z_sun<=met_arr[8])].f_gw.values)[::1], levels=levels,fill=False, 
               ax=ax[2,1], color=colors[1], zorder=3, linewidths=2.5)

    ax[2,2].scatter(y=utils.chirp_mass(COplot.mass_1.values*u.M_sun, 
                                  COplot.mass_2.values*u.M_sun),
               x=np.log10(COplot.f_gw.values), color='xkcd:light grey', zorder=0.)

    sb.kdeplot(y=utils.chirp_mass(CO.loc[(CO.met*Z_sun>=met_arr[-2])].mass_1.values*u.M_sun, 
                                  CO.loc[(CO.met*Z_sun>=met_arr[-2])].mass_2.values*u.M_sun)[::100],
               x=np.log10(CO.loc[(CO.met*Z_sun>=met_arr[-2])].f_gw.values)[::100], levels=levels,fill=False, 
               ax=ax[2,2], color=colors[3], zorder=3, linewidths=2.5)


    ax[3,0].scatter(y=utils.chirp_mass(ONeplot.mass_1.values*u.M_sun, 
                                  ONeplot.mass_2.values*u.M_sun),
               x=np.log10(ONeplot.f_gw.values), color='xkcd:light grey', zorder=0.)

    sb.kdeplot(y=utils.chirp_mass(ONe.loc[ONe.met*Z_sun<=met_arr[1]].mass_1.values*u.M_sun, 
                                  ONe.loc[ONe.met*Z_sun<=met_arr[1]].mass_2.values*u.M_sun)[::1],
               x=np.log10(ONe.loc[ONe.met*Z_sun<=met_arr[1]].f_gw.values)[::1], levels=levels,fill=False, 
               ax=ax[3,0], color=colors[0], zorder=3, linewidths=2.5)


    ax[3,1].scatter(y=utils.chirp_mass(ONeplot.mass_1.values*u.M_sun, 
                                  ONeplot.mass_2.values*u.M_sun),
               x=np.log10(ONeplot.f_gw.values), color='xkcd:light grey', zorder=0.)

    sb.kdeplot(y=utils.chirp_mass(ONe.loc[(ONe.met*Z_sun>=met_arr[7])&(ONe.met*Z_sun<=met_arr[8])].mass_1.values*u.M_sun, 
                                  ONe.loc[(ONe.met*Z_sun>=met_arr[7])&(ONe.met*Z_sun<=met_arr[8])].mass_2.values*u.M_sun)[::1],
               x=np.log10(ONe.loc[(ONe.met*Z_sun>=met_arr[7])&(ONe.met*Z_sun<=met_arr[8])].f_gw.values)[::1], levels=levels,fill=False, 
               ax=ax[3,1], color=colors[1], zorder=3, linewidths=2.5)

    ax[3,2].scatter(y=utils.chirp_mass(ONeplot.mass_1.values*u.M_sun, 
                                  ONeplot.mass_2.values*u.M_sun),
               x=np.log10(ONeplot.f_gw.values), color='xkcd:light grey', zorder=0.)

    sb.kdeplot(y=utils.chirp_mass(ONe.loc[(ONe.met*Z_sun>=met_arr[-2])].mass_1.values*u.M_sun, 
                                  ONe.loc[(ONe.met*Z_sun>=met_arr[-2])].mass_2.values*u.M_sun)[::10],
               x=np.log10(ONe.loc[(ONe.met*Z_sun>=met_arr[-2])].f_gw.values)[::10], levels=levels,fill=False, 
               ax=ax[3,2], color=colors[3], zorder=3, linewidths=2.5)

    from matplotlib.lines import Line2D
    custom_lines = [Line2D([0], [0], color='xkcd:light grey', lw=4),
                    Line2D([0], [0], color=colors[0], lw=4),
                    Line2D([0], [0], color=colors[1], lw=4),
                    Line2D([0], [0], color=colors[2], lw=4)]

    ax[0,0].legend([custom_lines[0], custom_lines[1]], ['All Z', 'Z=0.0001'], loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=4, borderaxespad=0, frameon=False, 
              fontsize=18)

    ax[0,1].legend([custom_lines[0], custom_lines[3]], ['All Z', 'Z={}'.format(met_arr[8])], loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=4, borderaxespad=0, frameon=False, 
              fontsize=18)

    ax[0,2].legend([custom_lines[0], custom_lines[2]], ['All Z', 'Z=03'], loc='lower left', bbox_to_anchor= (0.0, 1.01), ncol=4, borderaxespad=0, frameon=False, 
              fontsize=18)

    for i in range(4):
        ax[i,0].set_ylabel('Chirp Mass (M$_\odot$)', fontsize=20)
        ax[i,1].set_yticklabels('')
        ax[i,2].set_yticklabels('')
        #ax[i,1].set_yticks([])
        #ax[i,2].set_yticks([])

    for i in range(3):
        ax[3,i].set_xlabel(r'Log$_{10}$(f$_{\rm{GW}}$/Hz)', fontsize=20)
        ax[i,0].set_xticklabels('')
        ax[i,1].set_xticklabels('')
        ax[i,2].set_xticklabels('')
        ax[3,i].set_xticks([-4.25, -3.75, -3.25, -2.75])
        ax[3,i].set_xticklabels(['-4.25', '-3.75', '-3.25', '-2.75'])
        ax[0,i].set_ylim(0.175, 0.375)
        ax[2,i].set_ylim(0.3, 1.05)
        ax[1,i].set_ylim(0.2, 0.6)
        ax[3,i].set_ylim(0.3, 1.1)
        ax[0,i].text(0.85, 0.85, 'He + He', fontsize=18, horizontalalignment='center', 
                     transform=ax[0,i].transAxes)
        ax[1,i].text(0.85, 0.85, 'CO + He', fontsize=18, horizontalalignment='center', 
                     transform=ax[1,i].transAxes)
        ax[2,i].text(0.85, 0.85, 'CO + CO', fontsize=18, horizontalalignment='center', 
                     transform=ax[2,i].transAxes)
        ax[3,i].text(0.85, 0.85, 'ONe + X', fontsize=18, horizontalalignment='center', 
                     transform=ax[3,i].transAxes)

    for i in range(4):
        for j in range(3):
            ax[i,j].set_xlim(-4.25, -2.5)
            ax[i,j].axvline(-4.0, color='xkcd:grey', ls='--', lw=2., zorder=1.)

    plt.subplots_adjust(hspace=0.06, wspace=0.03)

    plt.show(block=False)
    
    return