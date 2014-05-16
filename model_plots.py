import sfh_tests_multi_proc as sfh_tests
import os
import numpy as np
import matplotlib.pylab as plt
from TPAGBparams import snap_src, research_path
import brewer2mpl
import ResolvedStellarPops as rsp
import fileIO
import galaxy_tests


color_scheme = ['#d73027', '#fc8d59', '#fee090', '#669966', '#e0f3f8', '#4575b4']
fontlarge = 24
fontmid = 20
fontsmall = 18
def translate_model_name(model, small=False):
    if 'oct' in model.lower():
        name = 'R75'
    if 'nov13eta' in model.lower():
        name = '\eta=0'
    if model.lower() == 'nov13':
        name = 'mSC05'
    if 'feb' in model.lower():
        name = 'FEB14'
    if small is True:
        new_model = r'$\dot M_{\rm pre\!-\!dust}^{\rm %s}$' % name
    else:
        new_model = r'$\dot{M}_{\rm pre\!-\!dust}=%s$' % name
    return new_model


def compare_agb_lifetimes():
    import glob
    track_loc = research_path + \
    'TP-AGBcalib/AGBTracks/plots_for_paperI/agbz001_3dup/'

    models = ['NOV13', 'NOV13eta0', 'OCT13']
    model_name = translate_model_name('nov13')

    # these two have to line up:
    search_formats = ['*dL0.0*', '*dL0.50*', '*dL2*']
    labels = [r'%s' % model_name,
              r'%s: $0.5\lambda$' % model_name,
              r'%s: $2\lambda$'  % model_name]

    track_sets = [rsp.fileIO.get_files(track_loc, sf) for sf in search_formats]
    cols1 =  ['k', '#d73027', '#fee090', '#e0f3f8', '#91bfdb', '#4575b4']
    bmap = brewer2mpl.get_map('Blues', 'Sequential', 9)
    cols2 = bmap.mpl_colors[3::2]

    fig, axs = plt.subplots(ncols=2, figsize=(10, 5), sharex=True, sharey=True)

    for i, track_set in enumerate(track_sets):
        tracks = np.array([fileIO.get_numeric_data(t) for t in track_set])
        tracks = tracks[np.argsort([t.mass for t in tracks])]
        masses = np.array([t.mass for t in tracks])
        taus = np.array([np.sum(t.data_array['dt']) for t in tracks])
        plt_kw = {'lw': 3, 'label': labels[i], 'color': cols1[i]}
        if i == 0:
            for ax in axs:
                ax.plot(masses, taus/1e6, lw=4, color='k')
                plt_kw['color'] = cols2[0]
                ax.plot(masses, taus/1e6, **plt_kw)
        else:
            axs[1].plot(masses, taus/1e6, lw=4, color='k')
            axs[1].plot(masses, taus/1e6, **plt_kw)

    for j in range(len(models)):
        model_name = models[j].replace('.dat', '').split('_')[-1]
        if models[j].lower() == 'nov13':
            continue

        base = research_path + \
            'TP-AGBcalib/AGBTracks/CAF09/S_%s' % model_name

        agb_track_loc = os.path.join(base, glob.glob1(base, '*0.001*')[0])

        track_names = [os.path.join(agb_track_loc, a)
                       for a in os.listdir(agb_track_loc)
                       if a.startswith('agb_') and not 'b_1.75' in a
                       and not 'b_1.80' in a]
        tracks = [fileIO.get_numeric_data(t) for t in track_names]
        tracks = [t for t in tracks if not t == -1 and t.data_array.size > 1]
        tracks = np.array(tracks)[np.argsort([t.mass for t in tracks])]

        masses = np.array([t.mass for t in tracks])
        taus = np.array([np.sum(t.data_array['dt']) for t in tracks])

        model_name = translate_model_name(model_name)
        plt_kw = {'lw': 3, 'label': model_name, 'color': cols2[j]}
    axs[0].plot(masses, taus/1e6, lw=4, color='k')
    axs[0].plot(masses, taus/1e6, **plt_kw)

    for ax in axs:
        ax.legend(loc=0, frameon=False, fontsize=fontsmall)
        ax.set_xlim(1, 2.95)
        ax.set_ylim(.25, 3.7)
        ax.set_xlabel(r'${\rm Initial\ Mass\ (M_\odot)}$', fontsize=fontlarge)
        ax.tick_params(labelsize=fontmid)

    fig.subplots_adjust(left=0.08, right=0.98, bottom=0.15, top=0.95,
                        wspace=0.01)
    axs[0].set_ylabel(r'${\rm Lifetime\ (Myr)}$', fontsize=fontlarge)

    plt.savefig('lambda_plot.png', dpi=150)
    return axs


def agb_lifetimes(models, z=0.002):
    import glob
    tauss = []
    btauss = []
    for j in range(len(models)):
        print models[j]
        print
        fig, ax = plt.subplots()
        fig2, ax2 = plt.subplots()
        model_name = models[j].replace('.dat', '').split('_')[-1]
        agb_track_loc = research_path + \
            'TP-AGBcalib/AGBTracks/CAF09/S_%s/' % model_name
        base = research_path + \
            'TP-AGBcalib/AGBTracks/CAF09/S_%s' % model_name
        if z == 'all':
            zs = np.array([d.split('_')[1].replace('Z','')
                           for d in glob.glob1(base, '*')], dtype=float)
        else:
            zs = [z]
        zs = np.array([i for i in zs if i <= 0.008])
        if len(zs) > 8.:
            zs = zs[::2]
        print zs
        cnum = np.max([len(zs), 3])
        bmap = brewer2mpl.get_map('Blues', 'Sequential', cnum + 1)
        cols = bmap.mpl_colors[1:]
        for i, z in enumerate(np.sort(zs)):
            try:
                agb_track_loc = os.path.join(base, glob.glob1(base, '*%g*' % z)[0])
            except IndexError:
                print 'no Z=%g tracks in %s' % (z, base)
                continue
            if not os.path.isdir(agb_track_loc) is True:
                print model_name, 'no agb tracks found'
            model_name = translate_model_name(models[j])
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
            with open('tpagb_lifetimes_S_%s_Z%g.dat' % (models[j], z), 'w') as out:
                out.write('# mass tpagb_tau tpagb_tau_bright \n')
                np.savetxt(out, np.array([masses, taus, btaus]).T, fmt='%.3f')


        for ax in [ax, ax2]:
            ax.set_xlabel('${\\rm Initial\ Mass\ (M_\odot)}$', fontsize=fontlarge)
            ax.set_ylabel('${\\rm Lifetime\ (Myr)}$', fontsize=fontlarge)
            ax.legend(loc=0, frameon=False)
            ax.set_xlim(0, 5)
            ax.set_ylim(0, 5)
        ax2.annotate('$\log L/L_\odot > 3.4$', (0.03, 0.97), xycoords='axes fraction',
                     fontsize=fontlarge, va='top')
        ax2.set_ylabel(ax2.get_ylabel().replace('(Myr)', '(Myr)'))
        #ax2.set_ylabel('${\\rm Pre\!-\!Flash\ Core\ Mass\ (M_\odot)}$', fontsize=24)

        fig.savefig('tpagb_lifetime_%s.png' % (models[j]), dpi=150)
        fig2.savefig('tpagb_lifetime_bright_%s.png' % (models[j]), dpi=150)

    return


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
    outfile_dir = snap_src + '/models/varysfh/match-hmc/'
    vSFH, vsfh_kws = sfh_tests.prepare_vsfh_run([target], ['cmd_input_CAF09_S_NOV13.dat'], 50,
                                                 vsfh_kw={'outfile_loc': outfile_dir,
                                                          'extra_str': ''},
                                                          default_kw=None)
    pl = sfh_tests.Plotting(vSFH[0])
    cols = color_scheme
    cols.append('#9966cc')
    del cols[2]
    del cols[3]
    kw = {'trilegal_output': trilegal_output,
          'narratio': False,
          'add_stage_lfs': 'all',
          'plot_data': False,
          'plot_models': False,
          'cols': cols,
          'stage_lf_kw': {'lw': 3, 'label': translate_model_name('nov13')}}
    ax1, ax2 = pl.compare_to_gal(target, **kw)
    lims = load_plot_limits()
    row = lims[np.nonzero(lims['target'] == target)[0]]
    for ax, band in zip([ax1, ax2], ['opt', 'ir']):
        ax.set_xlim(row['%s_cmdmin' % band], row['%s_cmdmax' % band])
        ax.set_ylim(row['%s_lfmin' % band], row['%s_lfmax' % band])


    base_dir = snap_src + '/models/varysfh/match-hmc'
    model = 'caf09_s_nov13'
    file_loc = os.path.join(base_dir, target, model, 'mc')
    best_tri_out, = rsp.fileIO.get_files(file_loc, '*opt_best*')
    filter1 = sfh_tests.get_filter1(target)
    sgal = rsp.Galaxies.simgalaxy(best_tri_out, filter1=filter1,
                                  filter2='F814W')

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


def compare_mass_loss(masses=1.0, z=0.001, sets=['NOV13', 'OCT13', 'NOV13eta0'],
              paola=False):
    '''
    made to plot a comparison between several mass prescriptions.
    Labels for the plot are set up stupidly, maybe in in_dict or labels arg...
    '''
    from matplotlib.ticker import NullFormatter

    teff_max = 3.5
    track_files = None
    if paola is True:
    # hack to use specific tracks from paola
        track_dir = research_path + '/TP-AGBcalib/AGBTracks/plots_for_paperI/'
        file_end = '_Mc0.00_dMc0.00_Tbd6.40_L0.00_dL0.00_C0.00_Nr3.00_rates0_KOPv_KMOLv.dat'
        if masses == 2.0:
            track_files = \
             [track_dir + 'agb_2.00_Z0.00100000_Mdot50_eta0.00' + file_end,
              track_dir + 'agb_2.00_Z0.00100000_Mdot49_eta0.40' + file_end,
              track_dir + 'agb_2.00_Z0.00100000_Mdot48_eta8.00' + file_end,
              track_dir + 'agb_2.00_Z0.00100000_Mdot50_eta0.40' + file_end]
            teff_max = 3.4
        if masses == 1.0:
            track_files = \
             [track_dir + 'agb_1.00_Z0.00100000_Mdot50_eta0.00' + file_end,
              track_dir + 'agb_1.00_Z0.00100000_Mdot49_eta0.40' + file_end,
              track_dir + 'agb_1.00_Z0.00100000_Mdot48_eta8.00' + file_end,
              track_dir + 'agb_1.00_Z0.00100000_Mdot50_eta0.40' + file_end]
            teff_max = 3.4

    labels = ['$\\dot{M}_{\\rm{pre-dust}}=0.0$','$\\rm{R75}$',
          '$\\rm{SC05}$', '$\\rm{mSC05}$']

    if track_files is not None:
        nrows = len(track_files)
    else:
        nrows = len(sets)
    fig, axs = plt.subplots(nrows=nrows, ncols=2, figsize=(8, 8))
    anorm = 1e6
    xlab0 = '\\rm{Age}\ (10^6\ \\rm{yr})'
    ylab0 = '\log\ \dot{M}\ (\\rm{M_\odot\ yr}^{-1})'
    ylab1 = '\log\ L\ (L_\odot)'
    xlab1 = '\log\ T_{\\rm eff}\ (\\rm{K})'

    agb_tracks_dir = research_path + 'TP-AGBcalib/AGBTracks/CAF09'

    if type(masses) is not list:
        masses = [masses]
    cols = ['k']

    for j, mass in enumerate(masses):
        if track_files is None:
            tnames = []
            labels = []
            for tset in sets:
                label = translate_model_name(tset)
                direc = os.path.join(agb_tracks_dir, 'S_' + tset)
                direc, = [os.path.join(direc, d)
                        for d in os.listdir(direc) if str(z) in d]
                tname = rsp.fileIO.get_files(direc, 'agb_%.2f*' % mass)[0]
                tnames.append(tname)
                labels.append('$%s$' % label)
            tracks = [fileIO.get_numeric_data(t) for t in tnames]
        else:
            tracks = [fileIO.get_numeric_data(t) for t in track_files]

        for i in range(len(tracks)):
            axs[i][0].plot(tracks[i].data_array['ageyr']/anorm,
                           tracks[i].data_array['dMdt'],
                           label='$M=%g\ M_\odot$' % mass, lw=1, color=cols[j])
            axs[i][0].plot(tracks[i].data_array['ageyr'][tracks[i].cstar]/anorm,
                           tracks[i].data_array['dMdt'][tracks[i].cstar],
                           lw=1, color='darkred')
            axs[i][1].plot(tracks[i].data_array['T_star'],
                           tracks[i].data_array['L_star'],
                           label='$M=%g\ M_\odot$' % mass, lw=1, color=cols[j])
            axs[i][1].plot(tracks[i].data_array['T_star'][tracks[i].cstar],
                           tracks[i].data_array['L_star'][tracks[i].cstar],
                           lw=1, color='darkred')
            axs[i][0].annotate(labels[i], (0.03, 0.96), fontsize=fontlarge,
                               xycoords='axes fraction', va='top')
    axs[-1, 0].set_xlabel('$%s$' % xlab0, fontsize=fontlarge)
    axs[-1, 1].set_xlabel('$%s$' % xlab1, fontsize=fontlarge)
    plt.annotate('$%s$' % ylab0, (0.03, 0.5), fontsize=fontlarge, va='center',
               xycoords='figure fraction', rotation='vertical')
    plt.annotate('$%s$' % ylab1, (0.95, 0.5), fontsize=fontlarge, va='center',
               xycoords='figure fraction', rotation='vertical')

    [ax.yaxis.tick_right() for ax in axs.flatten()[1::2]]
    [ax.xaxis.set_major_formatter(NullFormatter())
     for ax in axs.flatten()[:-2]]

    # mass loss
    [ax.set_ylim(-11.5, -4.5) for ax in axs[:, 0]]

    # log l
    [ax.set_ylim(2.81, 4.25) for ax in axs[:, 1]]

    # log te
    [ax.set_xlim(3.66, teff_max) for ax in axs[:, 1]]

    # age Myr
    [ax.set_xlim(0, 2.45) for ax in axs[:, 0]]

    # top left plot only
    if paola is False:
        [ax.legend(loc=4, fontsize=fontlarge, frameon=False)
     for ax in [axs.flatten()[0]]]

    fig.subplots_adjust(wspace=0.02, hspace=0.02)
    plt.savefig('compare_massloss_M%g_Z%g.png' % (masses[0], z), dpi=150)
    return axs


def tpagb_mass_histograms(chi2_location='draft_run', band='opt', dry_run=True,
                          model='nov13', model_src='default', force=False,
                          cumsum=True):
    '''
    plot a histogram of the scaled number of tpagb stars for the best fitting
    model in each model chi2 file in the chi2_location

    all args besides chi2_location are passed to tpagb_masses.

    the trilegal output files can be in a difference location with directories
    model_src/target/model/mc see tpagb_masses.

    plot colors are fixed at 6: errors for more than 6 targets.
    '''
    if chi2_location == 'draft_run':
        chi2_location = snap_src + '/models/varysfh/match-hmc/'

    chi2files = rsp.fileIO.get_files(chi2_location, '*%s_*chi2.dat' % model)
    # I do gaussian chi2 too, not just poisson...
    chi2files = [c for c in chi2files if not 'gauss' in c][::-1]

    # get the tpagb masses
    (masses, norm) = zip(*[tpagb_masses(c, band=band, dry_run=dry_run,
                                        model_src=model_src, force=force)
                    for c in chi2files])
    norm = np.array(norm)

    ts = [os.path.split(c)[1].split('_')[3] for c in chi2files]
    targets = galaxy_tests.ancients()
    tinds = [ts.index(t.lower()) for t in targets]
    targets = np.array(ts)[tinds]
    labels = ['$%s$' % t.upper().replace('-DEEP', '').replace('-', '\!-\!')
              for t in targets]


    # get the hist made nice and stacked, and then scale it down.
    hists, bins, pps = plt.hist(masses, stacked=True, align='left',
                                histtype='step', bins=50, visible=False)
    plt.close()

    # scaled histograms
    norm_hists = [hists[i] * norm[i] for i in range(len(hists))]

    # actual plot... norm scales the histogram.
    # set up plot
    cols = color_scheme

    if cumsum is True:
        fig, ax = plt.subplots()
        axs = [ax]
    else:
        fig, axs = plt.subplots(nrows=len(hists), sharex=True)
    # mask 0 values so there is a vertical line on the plot

    for i in range(len(hists)):
        norm_hists[i][norm_hists[i]==0] = 1e-5
        if cumsum is True:
            yplot = np.cumsum(norm_hists[i]) / np.sum(norm_hists[i])
        else:
            yplot = norm_hists[i]
        axs[i].plot(bins[:-1], yplot, linestyle='steps-pre', color='grey', lw=4)
        axs[i].plot(bins[:-1], yplot, linestyle='steps-pre', color=cols[i],
                    lw=2, label=labels[i], alpha=.9)

        #ax.plot(bins[:-1], np.sum(norm_hists, axis=0), linestyle='steps-pre',
        #             color='darkgrey', lw=3, label=r'$\rm{Total}$')
        axs[i].tick_params(labelsize=fontmid)


        #ax.set_yscale('log')
        #ax.set_ylim(3, 10**3)

    axs[-1].set_xlabel(r'$\rm{Mass\ M_\odot}$', fontsize=fontlarge)
    if cumsum is True:
        axs[0].legend(loc=0, frameon=False)
        axs[0].set_ylabel(r'$\rm{Cumulative\ Fraction\ of\ {TP\!-\!AGB}\ Stars}$', fontsize=fontlarge)
        axs[0].set_xlim(0.6, 3)
        fname = 'tpagb_mass_hist_%s_%s.png' % (band, model)
    else:
        [ax.set_xlim(0.6, 3) for ax in axs]
        plt.annotate(r'$\rm{\#\ of\ {TP\!-\!AGB}\ Stars}$', (0.03, 0.5),
                     fontsize=fontlarge, va='center',
               xycoords='figure fraction', rotation='vertical')
        fname = 'tpagb_mass_hist_%s_%s_nocumsum.png' % (band, model)
        [ax.locator_params(axis='y', nbins=2) for ax in axs]
        fig.subplots_adjust(hspace=0.03)
    plt.savefig(fname, dpi=150)
    return axs


def tpagb_masses(chi2file, band='opt', model_src='default', dry_run=False,
         mass=True, old=False, force=False):
    '''
    using the chi2file run trilegal with the best fit sfh and return the
    normalization and tp-agb masses (scaled simulation)
    '''
    if model_src == 'default':
        model_src = snap_src + '/models/varysfh/'


    components = os.path.split(chi2file)[1].split('_')
    model = '_'.join(components[:3])
    target = components[3]
    model_loc = os.path.join(model_src, target, model, 'mc')
    if force is False:
        # read chi2 file
        chi2_data = rsp.fileIO.readfile(chi2file)

        # best fitting chi2 run in band
        ibest_fit = np.argmin(chi2_data['%s_chi2' % band])

        # should work out to isfr == ibest_fit, but just in case:
        isfr = chi2_data['sfr'][ibest_fit]

        # associated input file for best run
        tri_inp, = rsp.fileIO.get_files(model_loc, '*%03d.dat' % isfr)
        tri_outp = tri_inp.replace('.dat', '_%s_best.dat' % band).replace('inp', 'outp')
    else:
        tri_inp, = rsp.fileIO.get_files(model_loc, '*best.dat')
        tri_outp = tri_inp.replace('inp', 'outp')
    rsp.fileIO.ensure_file(tri_inp)
    # run trilegal with best run

    cmd_input = 'cmd_input_%s.dat' % model.upper()
    rsp.TrilegalUtils.run_trilegal(cmd_input, tri_inp, tri_outp,
                                   dry_run=dry_run)

    # load trilegal run and do the normalization (I should really save that...)
    # all that's needed here is opt_norm and ir_norm.
    filter1 = sfh_tests.get_filter1(target)
    ags=sfh_tests.load_default_ancient_galaxies()

    files = sfh_tests.FileIO()
    files.read_trilegal_catalog(tri_outp, filter1=filter1)
    files.load_data_for_normalization(target=target, ags=ags)
    files.load_trilegal_data()
    sopt_rgb, sir_rgb, sopt_agb, sir_agb = \
        sfh_tests.rgb_agb_regions(files.sgal, files.opt_offset,
                                  files.opt_trgb, files.opt_trgb_err,
                                  ags, files.ir_offset,
                                  files.ir_trgb, files.ir_trgb_err,
                                  files.opt_mag, files.ir_mag)
    opt_norm, ir_norm, opt_rgb, ir_rgb, opt_agb, ir_agb = \
        sfh_tests.normalize_simulation(files.opt_mag, files.ir_mag,
                                       files.nopt_rgb, files.nir_rgb,
                                       sopt_rgb, sir_rgb, sopt_agb,
                                       sir_agb)
    print target, opt_norm, ir_norm
    with open(tri_inp.replace('.dat', '_norms.dat'), 'w') as out:
        out.write('# opt_norm ir_norm\n')
        out.write('%.6f %.6f \n' % (opt_norm, ir_norm))

    if band == 'opt':
        norm = opt_norm
    if band == 'ir':
        norm = ir_norm

    # load tp-agb stars
    files.sgal.all_stages('TPAGB')
    # it's crazy to believe there is a tp-agb star bluer than the color cut
    # but to be consistent with the rest of the method, will do a color cut.
    cut_inds = files.__getattribute__('%s_color_cut' % band)
    itpagb = list(set(files.sgal.itpagb) & set(cut_inds))

    # grab TP-AGB masses

    if mass is True:
        mass = files.sgal.data.get_col('m_ini')
        ret_val = mass[itpagb]
    else:
        met = files.sgal.data.get_col('[M/H]')
        if old is True:
            olds, = np.nonzero(files.sgal.data.get_col('logAge') > 8.5)
            itpagb = list(set(files.sgal.itpagb) & set(cut_inds) & set(olds))
        ret_val = met[itpagb]
    return ret_val, norm


def trilegal_metals(chi2_location='draft_run', band='opt', dry_run=False,
                    model='nov13', model_src='default', old=False, feh=False):
    if chi2_location == 'draft_run':
        chi2_location = snap_src + '/models/varysfh/match-hmc/'

    chi2files = rsp.fileIO.get_files(chi2_location, '*%s_*chi2.dat' % model)
    # I do gaussian chi2 too, not just poisson...
    chi2files = [c for c in chi2files if not 'gauss' in c][::-1]

    # get the tpagb masses
    (mhs, norm) = zip(*[tpagb_masses(c, band=band, dry_run=dry_run,
                                        model_src=model_src, mass=False,
                                        old=old) for c in chi2files])
    ts = [os.path.split(c)[1].split('_')[3] for c in chi2files]
    targets = galaxy_tests.ancients()
    tinds = [ts.index(t.lower()) for t in targets]
    targets = np.array(ts)[tinds]
    from ResolvedStellarPops.convertz import convertz
    if feh is True:
        ind = 4
    else:
        ind = 1
    zs = np.array([convertz(mh=i)[ind] for i in mhs])
    for i, target in enumerate(targets):
        print '%.4f %.4f %.4f %s ' % (np.min(zs[i]), np.median(zs[i]),
                                      np.max(zs[i]), target)


def plot_random_sfhs(targets='ancients'):
    '''
    plot the random sfr arrays with the data sfr arrays.
    Hard coded file locations. So read the code.
    '''
    targets = galaxy_tests.load_targets(targets)

    for target in targets:
        target = target.lower()
        # it doesn't matter which model so it is hard coded...
        sfr_file_loc = os.path.join(snap_src, 'models', 'varysfh', target,
                                    'caf09_s_nov13', 'mc')
        hmc_file_loc = os.path.join(snap_src, 'data', 'sfh_parsec')

        target = target.replace('-deep', '')
        outfile = '%s_random_sfr.png' % target
        hmc_file, = rsp.fileIO.get_files(hmc_file_loc, '%s*sfh' % target)

        sfh = sfh_tests.StarFormationHistories(hmc_file, 'match-hmc',
                                               sfr_file_loc=sfr_file_loc,
                                               sfr_file_search_fmt='*sfr')

        sfh.plot_sfh('sfr', plot_random_arrays_kw={'from_files': True},
                     outfile=outfile, zoom=True)
    return
