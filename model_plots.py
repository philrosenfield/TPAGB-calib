import sfh_tests
import os
import numpy as np
import matplotlib.pylab as plt
from TPAGBparams import snap_src, research_path
import brewer2mpl
import ResolvedStellarPops as rsp
import fileIO

def translate_model_name(model):
    if 'oct' in model.lower():
        name = '\\eta_R'
    if 'nov13eta' in model.lower():
        name = '\\eta=0'
    if model.lower() == 'nov13':
        name = 'SC05^\\prime'
    new_model = '$\\dot{M}_{%s}$' % name
    return new_model

def agb_lifetimes(models, z=0.002):
    import fileIO
    import glob
    tauss = []
    btauss = []
    fig, ax = plt.subplots()
    fig2, ax2 = plt.subplots()
    for i in range(len(models)):
        model_name = models[i].replace('.dat', '').split('_')[-1]
        agb_track_loc = research_path + 'TP-AGBcalib/AGBTracks/CAF09/S_%s/' % model_name
        base = research_path + 'TP-AGBcalib/AGBTracks/CAF09/S_%s' % model_name
        if z == 'all':
            zs = np.array([d.split('_')[1].replace('Z','') for d in glob.glob1(base, '*')], dtype=float)
        else:
            zs = [z]
        zs = np.array([i for i in zs if i <= 0.008])
        if len(zs) > 8.:
            zs = zs[::2]
        cnum = np.max([len(zs), 3])
        bmap = brewer2mpl.get_map('Blues', 'Sequential', cnum + 1)
        cols = bmap.mpl_colors[1:]
        for i, z in enumerate(np.sort(zs)):
            agb_track_loc = os.path.join(base, glob.glob1(base, '*%g*' % z)[0])

            if not os.path.isdir(agb_track_loc) is True:
                print model_name, 'no agb tracks found'
            model_name = translate_model_name(model_name)
            agb_track_names = [os.path.join(agb_track_loc, a)
                               for a in os.listdir(agb_track_loc)
                               if a.startswith('agb_')]
            tracks = [fileIO.get_numeric_data(agb_track)
                      for agb_track in agb_track_names]
            tracks = [t for t in tracks if not t == -1 and t.data_array.size > 1]
            masses = np.array([t.mass for t in tracks])
            sort = np.argsort(masses)
            masses = masses[sort]
            tracks = np.array(tracks)[sort]
            logls = np.array([t.get_col('L_star') for t in tracks])
            brights = np.array([np.nonzero(logl > 3.4)[0] for logl in logls])
            #m_cs = np.array([t.get_col('M_c')[0] for t in tracks])
            #ax2.plot(masses, m_cs, lw=2, color='black')

            taus = np.array([np.sum(t.data_array['dt']) for t in tracks])
            btaus = np.array([np.sum(t.data_array['dt'][b])
                              for t, b in zip(tracks, brights)])
            tauss.append(taus)
            btauss.append(btaus)
            plt_kw = {'lw': 3, 'label': '$Z=%g$' % z, 'color': cols[i]}
            ax.plot(masses, taus/1e6, lw=4, color='k')
            ax2.plot(masses, btaus/1e6, lw=4, color='k')

            ax.plot(masses, taus/1e6, **plt_kw)
            ax2.plot(masses, btaus/1e6, **plt_kw)
            with open('tpagb_lifetimes_S_%s_Z%g.dat' % (model_name, z), 'w') as out:
                out.write('# mass tpagb_tau tpagb_tau_bright \n')
                np.savetxt(out, np.array([masses, taus, btaus]).T, fmt='%.3f')

        #ax.fill_between(masses, tauss[0]/1e6, tauss[2]/1e6, alpha=0.1, color='grey')
        #ax2.fill_between(masses, btauss[0]/1e6, btauss[2]/1e6, alpha=0.1, color='grey')


        for ax in [ax, ax2]:
            ax.set_xlabel('${\\rm Initial\ Mass\ (M_\odot)}$', fontsize=20)
            ax.set_ylabel('${\\rm Lifetime\ (Myr)}$', fontsize=20)
            ax.legend(loc=0, frameon=False)
            ax.set_xlim(0, 5)
            ax.set_ylim(0, 5)
        ax2.set_ylabel(ax2.get_ylabel().replace('(Myr)', 'L>3.4L_\odot\ (Myr)'))
        #ax2.set_ylabel('${\\rm Pre\!-\!Flash\ Core\ Mass\ (M_\odot)}$', fontsize=20)
        fig.savefig('tpagb_lifetime_Z%g.png' % z, dpi=150)
        fig2.savefig('tpagb_lifetime_Z%g_bright.png' % z, dpi=150)


def load_plot_limits(filename='default'):
    if filename == 'default':
        filename = snap_src + '/tables/cmd_plot_limits.dat'
    dtype = [('target', '|S16'),
             ('opt_cmdmin', '<f8'),
             ('opt_cmdmax', '<f8'),
             ('opt_lfmin', '<f8'),
             ('opt_lfmax', '<f8'),
             ('ir_cmdmin', '<f8'),
             ('ir_cmdmax', '<f8'),
             ('ir_lfmin', '<f8'),
             ('ir_lfmax', '<f8'),
             ('opt_offset', '<f8'),
             ('ir_offset', '<f8')]
    lims = np.genfromtxt(filename, dtype=dtype)
    return lims.view(np.recarray)


def plot_lf_with_stages(target, trilegal_output):
    pl = sfh_tests.Plotting()
    cols = ['navy', 'darkgreen', 'darkgreen', 'darkgreen', 'purple',
            'darkred']
    kw = {'outfile_loc': '/home/phil/Dropbox/research/varysfh/',
          'trilegal_output': trilegal_output,
          'narratio': False,
          'hist_it_up': False,
          'add_stage_lfs': 'all',
          'plot_data': False,
          'cols': cols,
          'stage_lf_kw': {'lw': 3}}
    ax1, ax2 = pl.compare_to_gal(target, **kw)
    lims = load_plot_limits()['target'==target]
    for ax, band in zip([ax1, ax2], ['opt', 'ir']):
        ax.set_xlim(lims['%s_xmin' % band], lims['%s_xmax' % band])
        ax.set_ylim(lims['%s_ymin' % band], lims['%s_ymax' % band])

    return ax1, ax2


def plot_lfs():
    pl = sfh_tests.Plotting()
    outfile_loc='/home/phil/Dropbox/research/varysfh/'
    cmd_inputs = ['CAF09_S_NOV13'.lower(),
                  'CAF09_S_NOV13eta0'.lower(),
                  'CAF09_S_OCT13'.lower()]
    targets = ['ddo71', 'hs117', 'kkh37', 'ngc2976-deep', 'ngc404', 'ddo78']
    one_plot = True
    lims = load_plot_limits()
    for target in targets:
        print target
        if one_plot is True:
            fig, axs = plt.subplots(ncols=2, figsize=(12,6))
            plt.subplots_adjust(right=0.95, left=0.05, wspace=0.1)
            cols = ['black', 'navy', 'darkgreen']
            narratio=False
        else:
            cols = ['black'] * len(cmd_inputs)
            axs = None
            narratio=True
        plot_data = True
        for i, cmd_input in enumerate(cmd_inputs):
            if i > 0 and one_plot is True:
                plot_data = False
            narratio_file_name = os.path.join(outfile_loc,
                                              '%s_%s_narratio.dat' %
                                              (cmd_input, target.lower()))
            opt_lf_file = os.path.join(outfile_loc,
                                       '%s_%s_opt_lf.dat' %
                                       (cmd_input, target.lower()))
            ir_lf_file =  os.path.join(outfile_loc,
                                       '%s_%s_ir_lf.dat' %
                                       (cmd_input, target.lower()))

            ax1, ax2 = pl.compare_to_gal(target, opt_lf_file=opt_lf_file,
                                         ir_lf_file=ir_lf_file,
                                         hist_it_up=False, outfile_loc=outfile_loc,
                                         narratio_file_name=narratio_file_name,
                                         extra_str=cmd_input.split('_')[-1]+'_',
                                         axs=axs, plt_kw={'color': cols[i]},
                                         narratio=narratio, plot_data=plot_data)
            #ax1.set_title(cmd_input.replace('_', '\ '))
            lab = cmd_input.split('_')[-1]
            [ax.plot([0,0], [0,0], lw=3, color=cols[i], label='$%s$' % lab) for ax in [ax1, ax2]]
            row = lims[lims['target'] == target]

            ax1.set_xlim(row['opt_xmin'], row['opt_xmax'])
            ax2.set_xlim(row['ir_xmin'], row['ir_xmax'])
            ax1.set_ylim(row['opt_ymin'], row['opt_ymax'])
            ax2.set_ylim(row['ir_ymin'], row['ir_ymax'])
            [ax.legend(loc=0, frameon=False) for ax in [ax1, ax2]]
            figtitle = '%s%s_lfs.png' % (cmd_input.split('_')[-1]+'_', target)
            outfile = os.path.join(outfile_loc, figtitle)


def compare_mass_loss(mass=1.0, z=0.002, sets=['NOV13', 'OCT13', 'NOV13eta0']):
    '''
    made to plot a comparison between several mass prescriptions.
    Labels for the plot are set up stupidly, maybe in in_dict or labels arg...
    '''
    fig, axs = plt.subplots(nrows=len(sets), ncols=2, figsize=(10, 8))

    anorm = 1e6
    xlab0 = 'Age\ (10^6\ \\rm{yr})'
    ylab0 = '\log\ \dot{M}\ (\\rm{M_\odot\ yr}^{-1})'
    ylab1 = '\log\ L\ (L_\odot)'
    xlab1 = '\log\ T_{\\rm eff}\ (\\rm{K})'

    agb_tracks_dir = research_path + 'TP-AGBcalib/AGBTracks/CAF09'
    direcs = []
    labels = []
    for tset in sets:
        label = translate_model_name(tset)
        direc = os.path.join(agb_tracks_dir, 'S_' + tset)
        direc, = [os.path.join(direc, d)
                for d in os.listdir(direc) if str(z) in d]
        direc = rsp.fileIO.get_files(direc, 'agb_%.2f*' % mass)[0]
        direcs.append(direc)
        labels.append('$%s$' % label)
    tracks = [fileIO.get_numeric_data(t) for t in direcs]

    for i in range(len(tracks)):
        axs[i][0].plot(tracks[i].data_array['ageyr']/anorm,
                       tracks[i].data_array['dMdt'],
                       label=labels[i], lw=1, color='k')
        axs[i][1].plot(tracks[i].data_array['T_star'],
                       tracks[i].data_array['L_star'],
                       label=labels[i], lw=1, color='k')

    axs[-1, 0].set_xlabel('$%s$' % xlab0, fontsize=20)
    axs[-1, 1].set_xlabel('$%s$' % xlab1, fontsize=20)
    axs[1, 0].set_ylabel('$%s$' % ylab0, fontsize=20)
    axs[1, 1].set_ylabel('$%s$' % ylab1, fontsize=20)

    #ax.set_ylabel('$%s$' % ylab, fontsize=20)
    #ax.text(.03, .05, '$Z=%g\ M=%g\ M_\odot$' % (z, mass), fontsize=20,
    #        transform=ax.transAxes)
    [ax.tick_params(labelsize=16) for ax in axs.flatten()]
    [ax.set_ylim(-11.5, -4.5) for ax in axs[:, 0]]
    [ax.set_ylim(2.7, 4.2) for ax in axs[:, 1]]
    [ax.set_xlim(3.7, 3.48) for ax in axs[:, 1]]
    [ax.set_xlim(0, 2.5) for ax in axs[:, 0]]
    [ax.legend(loc=4, frameon=False, fontsize=16) for ax in axs.flatten()]
    fig.subplots_adjust(wspace=0.3)
    plt.savefig('compare_massloss.png', dpi=150)
    return ax
