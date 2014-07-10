import os
import matplotlib.pylab as plt
import numpy as np
import ResolvedStellarPops as rsp
from sfhs import star_formation_histories
__all__ = ['Plotting']


def plot_random_sfhs(vsfh):
    '''
    plot the random sfr arrays with the data sfr arrays.
    Hard coded file locations. So read the code.
    '''
    sfr_file_loc = vsfh.outfile_loc
    hmc_file_loc = vsfh.base

    outfile = os.path.join(sfr_file_loc, '%s_random_sfr.png' % vsfh.target)
    hmc_file, = rsp.fileIO.get_files(hmc_file_loc, '%s*mcmc.zc' % vsfh.target)

    sfh = star_formation_histories(hmc_file, 'match-hmc',
                                 sfr_file_loc=sfr_file_loc,
                                 sfr_file_search_fmt='*sfr')

    sfh.plot_sfh('sfr', plot_random_arrays_kw={'from_files': True},
                 outfile=outfile)
    return

class Plotting(object):
    def __init__(self, vsfh):
        key = ['target', 'agb_mod', 'fnames', 'offsets', 'trgb_excludes',
               'nopt_rgb', 'nir_rgb', 'nopt_agb', 'nir_agb', 'opt_trgb',
               'ir_trgb']
        [self.__setattr__(k, vsfh.__getattribute__(k)) for k in key]
        self.opt_lf_file, = [f for f in self.fnames if 'opt_lf' in f]
        self.ir_lf_file, = [f for f in self.fnames if 'ir_lf' in f]
        self.narratio_file, = [f for f in self.fnames if 'narratio' in f]
        self.vsfh = vsfh

    def plot_lf_file(self, opt_lf_file, ir_lf_file, axs=None, plt_kw=None,
                     opt_limit=None, ir_limit=None):
        '''needs work, but: plot the lf files.'''
        # set up the plot
        plt_kw = plt_kw or {}
        plt_kw = dict({'linestyle': 'steps-mid', 'color': 'black',
                       'alpha': 0.2}.items() + plt_kw.items())
        label = '%s' % self.agb_mod.split('_')[-1]
        plt_kw_lab = dict(plt_kw.items() + {'label': label}.items())
        if axs is None:
            fig, (axs) = plt.subplots(ncols=2, figsize=(12, 6))
            plt.subplots_adjust(right=0.95, left=0.05, wspace=0.1)

        # these have like 50 histograms each
        opt_hists, opt_binss, opt_norms = self.vsfh.load_lf_file(self.opt_lf_file)
        ir_hists, ir_binss, ir_norms = self.vsfh.load_lf_file(self.ir_lf_file)

        for i, (hists, binss, limit, norms) in \
            enumerate(zip([opt_hists, ir_hists], [opt_binss, ir_binss],
                          [opt_limit, ir_limit], [opt_norms, ir_norms])):

            for j, (hist, bins, norm) in enumerate(zip(hists, binss, norms)):
                if j != 0:
                    kw = plt_kw
                else:
                    kw = plt_kw_lab
                if limit is not None:
                    inds, = np.nonzero(bins <= limit)
                    axs[i].plot(bins[inds], hist[inds] * norm, **kw)
                else:
                    axs[i].plot(bins, hist * norm, **kw)

        return axs

    def count_stars_from_hist(self, opt_hist, opt_bins, ir_hist, ir_bins):
        ratio_data = {}
        for i, (hist, bins, band) in enumerate(zip([opt_hist, ir_hist],
                                                   [opt_bins, ir_bins],
                                                   ['opt', 'ir'])):
            trgb = self.__getattribute__('%s_trgb' % band)
            trgb_err = self.__getattribute__('%s_trgb_err' % band)
            norm = self.__getattribute__('%s_offset' % band)
            irgb = rsp.math_utils.between(bins, norm,
                                          trgb + trgb_err * self.ags.factor[i])
            iagb = rsp.math_utils.between(bins,
                                          trgb - trgb_err * self.ags.factor[i],
                                          10.)

            nrgb = np.sum(hist[irgb])
            nagb = np.sum(hist[iagb])
            ratio_data['%s_ar_ratio' % band] = nagb / nrgb
            ratio_data['%s_ar_ratio_err' % band] = \
                rsp.utils.count_uncert_ratio(nagb, nrgb)
            ratio_data['n%s_rgb' % band] = nrgb
            ratio_data['n%s_agb'% band] = nagb
        return ratio_data

    def add_narratio_to_plot(self, ax, band, ratio_data):
        stext_kw = dict({'color': 'black', 'fontsize': 14, 'ha': 'center'}.items() +
                        rsp.graphics.GraphicsUtils.ann_kwargs.items())
        dtext_kw = dict(stext_kw.items() + {'color': 'darkred'}.items())

        nrgb = ratio_data['n%s_rgb' % band]
        nagb = ratio_data['n%s_agb'% band]
        dratio = nagb / nrgb
        dratio_err = rsp.utils.count_uncert_ratio(nagb, nrgb)

        #yval = 1.2  # text yloc found by eye, depends on fontsize
        stext_kw['transform'] = ax.transAxes
        dtext_kw['transform'] = ax.transAxes
        yval = 0.95
        xagb_val = 0.17
        xrgb_val = 0.5
        xratio_val = 0.83
        xvals = [xagb_val, xrgb_val, xratio_val]

        # simulated nrgb and nagb are the mean values
        srgb_text = '$\langle N_{\\rm RGB}\\rangle =%i$' % \
                    np.mean(ratio_data['n%s_rgb' % band])
        sagb_text = '$\langle N_{\\rm TP-AGB}\\rangle=%i$' % \
                    np.mean(ratio_data['n%s_agb' % band])

        # one could argue taking the mean isn't the best idea for
        # the ratio errors.
        sratio_text = '$f=%.3f\pm%.3f$' % \
                      (np.mean(ratio_data['%s_ar_ratio' % band]),
                       np.mean(ratio_data['%s_ar_ratio_err' % band]))

        drgb_text = '$N_{\\rm RGB}=%i$' % nrgb
        dagb_text = '$N_{\\rm TP-AGB}=%i$' % nagb
        dratio_text =  '$f = %.3f\pm%.3f$' % (dratio, dratio_err)

        textss = [[sagb_text, srgb_text, sratio_text],
                 [dagb_text, drgb_text, dratio_text]]
        kws = [stext_kw, dtext_kw]

        for kw, texts in zip(kws, textss):
            for xval, text in zip(xvals, texts):
                ax.text(xval, yval, text, **kw)
            yval -= .05  # stack the text
        return ax

    def plot_by_stage(self, ax1, ax2, add_stage_lfs='default', stage_lf_kw=None,
                      cols=None, trilegal_output=None, hist_it_up=False,
                      narratio=True):

        if add_stage_lfs == 'all':
            add_stage_lfs = ['MS', 'RGB', 'HEB', 'RHEB',
                             'BHEB', 'EAGB', 'TPAGB']
        if add_stage_lfs == 'default':
            add_stage_lfs = ['RGB', 'EAGB', 'TPAGB']

        nstages = len(add_stage_lfs)
        stage_lf_kw = stage_lf_kw or {}
        stage_lf_kw = dict({'linestyle': 'steps', 'lw': 2}.items() +
                            stage_lf_kw.items())
        if hasattr(stage_lf_kw, 'label'):
            stage_lf_kw['olabel'] = stage_lf_kw['label']
        if cols is None:
            cols = color_scheme

        # load the trilegal catalog if it is given, if it is given,
        # no LF scaling... need to save this info better. Currently only
        # in log files.
        if trilegal_output is not None:
            self.files.read_trilegal_catalog(trilegal_output,
                                             filter1=get_filter1(self.target))
            self.files.load_trilegal_data()
            self.opt_norm = 1.
            self.ir_norm = 1.

        self.files.load_data_for_normalization(target=self.target, ags=self.ags)
        assert hasattr(self.files, 'opt_mag'), \
            'Need opt_mag or trilegal_output'

        for ax, mag, norm, sinds, bins in \
            zip([ax1, ax2],
                [self.files.sgal.data.get_col('F814W')-self.files.opt_moffset,
                 self.files.sgal.data.get_col('F160W')-self.files.ir_moffset],
                [self.opt_norm, self.ir_norm],
                [self.files.opt_color_cut, self.files.ir_color_cut],
                [self.files.opt_bins, self.files.ir_bins]):

            #self.files.sgal.make_lf(mag, stages=add_stage_lfs, bins=bins,
            #                        inds=sinds, hist_it_up=hist_it_up)
            self.files.sgal.make_lf(mag, stages=add_stage_lfs, bins=bins,
                                    hist_it_up=hist_it_up)
            #import pdb
            #pdb.set_trace()
            k = 0
            for i in range(nstages):
                istage = add_stage_lfs[i].lower()
                try:
                    hist = self.files.sgal.__getattribute__('i%s_lfhist' %
                                                            istage)
                except AttributeError:
                    continue
                # combine all HeB stages into one for a cleaner plot.
                if istage == 'heb':
                    hist = \
                    np.sum([self.files.sgal.__getattribute__('i%s_lfhist' %
                                                             jstage.lower())
                            for jstage in add_stage_lfs
                            if 'heb' in jstage.lower()], axis=0)
                elif 'heb' in istage:
                    continue

                bins = \
                self.files.sgal.__getattribute__('i%s_lfbins' %
                                                 istage)
                stage_lf_kw['color'] = cols[k]
                k += 1
                stage_lf_kw['label'] = '$%s$' % istage.upper().replace('PA', 'P\!-\!A')
                norm_hist = hist # * norm
                norm_hist[norm_hist < 3] = 3
                ax.plot(bins[:-1], norm_hist, **stage_lf_kw)

        sopt_hist, sopt_bins = self.files.sgal.make_lf(self.files.opt_mag,
                                                      bins=self.files.opt_bins,
                                                      hist_it_up=hist_it_up)
        sir_hist, sir_bins = self.files.sgal.make_lf(self.files.ir_mag,
                                                     bins=self.files.ir_bins,
                                                     hist_it_up=hist_it_up)
        # 3 is the plot limit...
        sopt_hist[sopt_hist < 3] = 3
        sir_hist[sir_hist < 3] = 3
        stage_lf_kw['color'] = 'black'

        if narratio is False:
            if not hasattr(stage_lf_kw, 'olabel'):
                lab = '$Total$'
                #if hasattr(self, 'agb_mod'):
                #    model = self.agb_mod.split('_')[-1]
                #    lab = model_plots.translate_model_name(model, small=True)
                stage_lf_kw['label'] = lab
            else:
                stage_lf_kw['label'] = stage_lf_kw['olabel']
        ax1.plot(sopt_bins[:-1], sopt_hist, **stage_lf_kw)
        ax2.plot(sir_bins[:-1], sir_hist, **stage_lf_kw)
        return ax1, ax2

    def plot_gal(self, opt_gal, ir_gal, ax1, ax2):
        dplot_kw = {'drawstyle': 'steps-mid', 'color': 'darkred', 'lw': 2,
                    'label': '$%s$' % self.target}
        # HACK to mask low values
        #opt_gal.hist[opt_gal.hist < 0.1] = 0.1
        #ir_gal.hist[ir_gal.hist < 0.1] = 0.1
        opt_gal.hist, opt_gal.bins = np.histogram(opt_gal.mag2,
                                                  bins=self.vsfh.opt_bins)
        ir_gal.hist, ir_gal.bins = np.histogram(ir_gal.mag2,
                                                bins=self.vsfh.ir_bins)
        opt_err = np.sqrt(opt_gal.hist)
        ir_err = np.sqrt(ir_gal.hist)
        ax1.errorbar(opt_gal.bins[1:], opt_gal.hist, yerr=opt_err,
                     **dplot_kw)
        ax2.errorbar(ir_gal.bins[1:], ir_gal.hist, yerr=ir_err,
                     **dplot_kw)
        return ax1, ax2

    def compare_to_gal(self, opt_gal, ir_gal, opt_limit, ir_limit,
                       narratio=True, no_agb=False, xlim=None, ylim=None,
                       extra_str='', cols=None, stage_lf_kw=None, axs=None,
                       plt_kw=None):
        '''
        Plot the LFs and galaxy LF.

        ARGS:
        narratio: overlay NRGB, NAGB, and NAGB/NRGB +/- err
        no_agb: plot the LF without AGB stars

        RETURNS:
        ax1, ax2: axes instances created for the plot.

        '''
        import pdb; pdb.set_trace()
        # plot lfs from simulations (and initialize figure)
        plt_kw = plt_kw or {}
        (ax1, ax2) = \
            self.plot_lf_file(self.opt_lf_file, self.ir_lf_file,
                              opt_limit=opt_limit,
                              ir_limit=ir_limit, axs=axs,
                              plt_kw=plt_kw)
        # plot galaxy data
        ax1, ax2 = self.plot_gal(opt_gal, ir_gal, ax1, ax2)

        # initialize add numbers to the plot
        if narratio is True:
            # count stars from the saved file
            ratio_data = rsp.fileio.readfile(self.narratio_file,
                                             string_column=0)
            # get the number ratios for the annotations

            mean_ratio = {}
            for key in ratio_data.dtype.names:
                if key == 'target':
                    continue
                mean_ratio[key] = np.mean(ratio_data[:][key])

        for i, (ax, gal, trgb) in enumerate(zip([ax1, ax2], [opt_gal, ir_gal],
                                                [self.opt_trgb, self.ir_trgb])):
            ax.set_yscale('log')
            if ylim is not None:
                ax.set_ylim(ylim)
            if xlim is not None:
                ax.set_xlim(xlim)

            yarr = np.linspace(*ax.get_ylim())
            # vertical lines around the trgb exclude region
            ax.fill_betweenx(yarr, trgb - self.trgb_excludes[i],
                             trgb + self.trgb_excludes[i],
                             color='black', alpha=0.1)

            ax.vlines(trgb, *ax.get_ylim(), color='black',
                      linestyle='--')

            loc = 4
            if narratio is False:
                loc = 0
            ax.legend(loc=loc)
            ax.set_xlabel('$%s$' % gal.filter2, fontsize=20)

            if narratio is True:
                # need to load the data nrgb and nagb, calculate the ratio
                # and error.
                #self.add_narratio_to_plot(ax, band, mean_ratio)
                pass

        ax1.set_ylabel(r'${\rm Number\ of\ Stars}$', fontsize=20)
        plt.tick_params(labelsize=16)
        outfile = '%s%s_lfs.png' % (self.opt_lf_file.split('opt_lf')[0][:-1],
                                    extra_str)
        plt.savefig(outfile, dpi=150)
        print 'wrote %s' % outfile
        return ax1, ax2

class DiagnosticPlots(Plotting):
    def __init__(self, vsfh):
        self.Plotting.__init__(vsfh)

    def plot_mass_met_table(self, opt_mass_met_file, ir_mass_met_file,
                            extra_str=''):
        from mpl_toolkits.axes_grid1 import ImageGrid
        from astroML.stats import binned_statistic_2d

        fig = plt.figure(figsize=(8, 8))
        grid = ImageGrid(fig, 111,
                         nrows_ncols=(2, 2),
                         axes_pad=.5,
                         add_all=True,
                         label_mode="all",
                         cbar_location="top",
                         cbar_mode="each",
                         cbar_size="7%",
                         cbar_pad="2%",
                         aspect=0)
        cmaps = [plt.cm.get_cmap('jet', 9), plt.cm.gray_r]
        #cmap =
        #cmap.set_bad('w', 1.)
        #fig, (axs) = plt.subplots(ncols=2, figsize=(8, 8), sharey=True)
        types = ['mean', 'count']
        k =-1
        for j in range(len(types)):
            for i, mass_met in enumerate([opt_mass_met_file,
                                          ir_mass_met_file]):
                k += 1
                with open(mass_met, 'r') as mmf:
                    lines = [l.strip() for l in mmf.readlines()
                             if not l.startswith('#')]

                mag = np.concatenate([np.array(l.split(), dtype=float)
                                      for l in lines[0::3]])
                mass = np.concatenate([np.array(l.split(), dtype=float)
                                       for l in lines[1::3]])
                mh = np.concatenate([np.array(l.split(), dtype=float)
                                     for l in lines[2::3]])

                N, xedges, yedges = binned_statistic_2d(mag, mass, mh,
                                                        types[j], bins=50)
                im = grid[k].imshow(N.T, origin='lower',
                               extent=[xedges[0], xedges[-1], yedges[0],
                                       yedges[-1]],
                               aspect='auto', interpolation='nearest',
                               cmap=cmaps[j])
                grid[k].cax.colorbar(im)
                #grid[i].cax.set_label('$[M/H]$')

        grid.axes_all[0].set_ylabel('${\\rm Mass}\ (M_\odot)$', fontsize=20)
        grid.axes_all[2].set_ylabel('${\\rm Mass}\ (M_\odot)$', fontsize=20)
        grid.axes_all[2].set_xlabel('$F814W$', fontsize=20)
        grid.axes_all[3].set_xlabel('$F160W$', fontsize=20)
        target = '_'.join(os.path.split(opt_mass_met_file)[1].split('_')[0:4])
        fig.suptitle('$%s$' % target.replace('_', '\ '), fontsize=20)
        plt.savefig('%s_mass_met%s.png' % (target, extra_str), dpi=150)
        return grid


def tpagb_mass_histograms(chi2_location='draft_run', band='opt', dry_run=True,
                         model='nov13', model_src='default'):
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
                                        model_src=model_src)
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

    fig, ax = plt.subplots()

    # mask 0 values so there is a vertical line on the plot
    for i in range(len(hists)):
        norm_hists[i][norm_hists[i]==0] = 1e-5
        yplot = np.cumsum(norm_hists[i]) / np.sum(norm_hists[i])
        #yplot = norm_hists[i]
        ax.plot(bins[:-1], yplot, linestyle='steps-pre', color='grey', lw=4)
        ax.plot(bins[:-1], yplot, linestyle='steps-pre', color=cols[i],
                lw=2, label=labels[i], alpha=.9)

    #ax.plot(bins[:-1], np.sum(norm_hists, axis=0), linestyle='steps-pre',
    #             color='darkgrey', lw=3, label=r'$\rm{Total}$')
    ax.legend(loc=0, frameon=False)

    #ax.set_yscale('log')
    #ax.set_ylim(3, 10**3)
    ax.set_xlim(0.6, 3)
    ax.set_xlabel(r'$\rm{Mass\ M_\odot}$', fontsize=20)
    ax.set_ylabel(r'$\rm{Cumulative\ Fraction\ of\ {TP\!-\!AGB}\ Stars}$', fontsize=20)
    plt.tick_params(labelsize=16)
    plt.savefig('tpagb_mass_hist_%s_%s.png' % (band, model), dpi=150)
    return ax


def tpagb_masses(chi2file, band='opt', model_src='default', dry_run=False,
		 mass=True, old=False):
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

    # read chi2 file
    chi2_data = rsp.fileIO.readfile(chi2file)

    # best fitting chi2 run in band
    ibest_fit = np.argmin(chi2_data['%s_chi2' % band])

    # should work out to isfr == ibest_fit, but just in case:
    isfr = chi2_data['sfr'][ibest_fit]

    # associated input file for best run
    tri_inp, = rsp.fileIO.get_files(model_loc, '*%03d.dat' % isfr)

    # run trilegal with best run
    tri_outp = tri_inp.replace('.dat', '_%s_best.dat' % band).replace('inp', 'outp')
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
