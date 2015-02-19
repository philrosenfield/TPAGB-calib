import ResolvedStellarPops as rsp
import pyfits
import numpy as np
import matplotlib.pylab as plt
#from pop_synth.stellar_pops import normalize_simulation, rgb_agb_regions
plt_labelsize = 20
tck_labelsize = 16


def rgb_agb_regions(offset, trgb_exclude, trgb, mag):
    # define RGB regions
    if offset >= 2.:
        low = offset
    else:
        low = offset + trgb
    mid = trgb + trgb_excludew

    # Recovered stars in simulated RGB region.
    srgb = rsp.starpop.stars_in_region(mag, low, mid)

    # define AGB regions
    mmid = trgb - trgb_exclude
    high = 10.

    # Recovered stars in simulated AGB region.
    sagb = rsp.starpop.stars_in_region(mag, mmid, high)
    return srgb, sagb


def normalize_simulation(mag, ndata_rgb, srgb, sagb):
    norm = ndata_rgb / float(len(srgb))
    #print 'Normalization: %f' % norm

    # random sample the data distribution
    rands = np.random.random(len(mag))
    ind, = np.nonzero(rands < norm)

    # scaled rgb: norm + in rgb
    rgb = list(set(ind) & set(srgb))

    # scaled agb
    agb = list(set(ind) & set(sagb))
    return norm, rgb, agb


def load_data(fname, field=None):
    gal = rsp.StarPop()
    gal.data = pyfits.getdata(fname)

    if hasattr(gal.data, 'F160W_GST'):
        gst, = np.nonzero((gal.data.F160W_GST == 'T') & (gal.data.F110W_GST == 'T') &
                          (gal.data.F475W_GST == 'T') & (gal.data.F814W_GST == 'T'))
        gal.gst = gst
    else:
        gal.gst = np.arange(len(gal.data.field('field')))

    if field is not None:
        finds, = np.nonzero(gal.data.field('field') == field)
        # If in overlapping fields, I added the field flag together.
        if field == 6 or field == 12:
            finds2, = np.nonzero(gal.data.field('field') == 18)
            finds = np.concatenate([finds,finds2])
        gal.finds = finds
    else:
        gal.finds = np.arange(len(gal.data.field('field')))
    return gal


def load_phot(StarPop, band, filter1=None, filter2=None, color_min=None):

    try:
        StarPop.icut = list(set(StarPop.gst) & set(StarPop.finds))
    except AttributeError:
        StarPop.icut = np.arange(len(StarPop.data[filter1]))
    mag_attr = '%s_mag' % band
    mag1_attr = '%s_mag1' % band
    color_attr = '%s_color' % band
    color_cut_attr = '%s_color_cut' % band
    if not None in [filter1, filter2]:
        color = StarPop.data[filter1] - StarPop.data[filter2]
        StarPop.__setattr__(color_attr, color)
        if color_min is not None:
            color_cut, = np.nonzero((color) > color_min)
            StarPop.cut = list(set(StarPop.icut) & set(color_cut))
            StarPop.__setattr__(color_cut_attr, color_cut)
        else:
            StarPop.cut = StarPop.icut
        StarPop.__setattr__(mag1_attr, StarPop.data[filter1][StarPop.cut])
        StarPop.__setattr__(mag_attr, StarPop.data[filter2][StarPop.cut])

    return StarPop


def load_model(fname):
    sgal = rsp.StarPop()

    # all the columns
    with open(fname, 'r') as f:
        col_keys = f.readline().strip().split()

    # the columns I want
    cols =   ['logAge',
              '[M/H]',
              'm_ini',
              'logL',
              'logTe',
              'mbol',
              'F475W',
              'F814W',
              'F110W',
              'F160W',
              'IRAC_3.6',
              'IRAC_4.5',
              'Mcore',
              'C/O',
              'Per',
              'mode',
              'logML',
              'Mact',
              'stage']
    usecols = [col_keys.index(c) for c in cols]
    sgal.data = np.genfromtxt(fname, usecols=usecols,
                              names=cols)

    return sgal


def make_lf(ghist, gbins, shist, sbins, gerrs=None, serrs=None,
            colors=['black', 'darkred'], label='', ax=None):

    if ax is None:
        fig, ax = plt.subplots(figsize=(8, 8))
    if not None in [gerrs, serrs]:
        plt_kw = {'drawstyle': 'steps-mid', 'color': colors[0],
                  'lw': 2}
        ax.errorbar(gbins[1:], ghist, yerr=gerrs, **plt_kw)
        plt_kw['color'] = colors[1]
        plt_kw['label'] = label
        plt_kw['alpha'] = 0.3
        ax.errorbar(sbins[1:], shist, yerr=serrs, **plt_kw)
    else:
        plt_kw = {'linestyle': 'steps-mid', 'color': colors[0],
                  'lw': 2}
        ax.plot(gbins[1:], ghist, **plt_kw)
        plt_kw['color'] = colors[1]
        plt_kw['label'] = label
        plt_kw['alpha'] = 0.3
        ax.plot(sbins[1:], shist, **plt_kw)

    ax.set_yscale('log')
    ax.tick_params(labelsize=16)

    ax.set_ylabel(r'${\rm Number\ of\ Stars}$', fontsize=20)
    return ax


def add_vlines_lf(ax, offset, trgb_exclude, trgb):
    yarr = np.linspace(*ax.get_ylim())
    ax.fill_betweenx(yarr, trgb - trgb_exclude, trgb + trgb_exclude,
                     color='black', alpha=0.1)
    ax.vlines(trgb, *ax.get_ylim(), color='black', linestyle='--')
    ax.vlines(offset, *ax.get_ylim(), color='black', linestyle='--')
    return ax

#offsets = [22.5, 19.2]
nir_offset = 19.2
opt_offset = 22.5
#trgb_excludes = [0.1, 0.15]
nir_trgb_exclude = 0.15
opt_trgb_exclude = 0.1

opt_trgb = 20.5
nir_trgb = 18.3
bands = ['opt', 'nir', 'ir']
color_mins = [3.5, 0.98, None]
# load data and model

base = '/home/rosenfield/research/TP-AGBcalib/PHAT/'

#gal_file = base + 'data/Spitzer/b21-6filt-cut-shallow-fields.fits'
#gfilters = ['F475W_VEGA', 'F814W_VEGA', 'F110W_VEGA', 'F160W_VEGA', None, None]
#sfilt1ers= ['F475W', 'F814W', 'F110W', 'F160W', None, None]

#gal_file = base + 'data/Spitzer/B21_phat_irac_fields.fits'
gal_file = base + 'data/Spitzer/B21_phat_irac_v2_fields.fits'
gfilters = [['F475W_VEGA', 'F814W_VEGA'], ['F110W_VEGA', 'F160W_VEGA'], ['IRAC1', 'IRAC2']]
sfilters = [['F475W', 'F814W'], ['F110W', 'F160W'], ['IRAC_36', 'IRAC_45']]

bins = np.arange(10, 23, 0.2)
fields = [6, 12, 15]
# the best sfh with aringer, the best sfh with loidl
extras = [''] #['.aringer', '']
colors = ['navy', 'darkred', 'orange']
agb_mods = ['nov13_pagb', 'mar13_pagb', 'oct13_pagb']
sgal_fmt = base + 'vary_sfh/output_M31-B21_3x6-0%02i_%s%s.dat'

for extra in extras:
    for field in fields:
        fig, axs = plt.subplots(ncols=3, figsize=(18, 8))
        for k, agb_mod in enumerate(agb_mods):
            # read the data files
            gal = load_data(gal_file, field=field)
            sgal_name = sgal_fmt % (field, agb_mod, extra)
            sgal = load_model(sgal_name)
            for i, filters in enumerate(gfilters):
                # load the photometery in each filter set
                gal = load_phot(gal, bands[i], filter1=filters[0],
                                filter2=filters[1], color_min=color_mins[i])

                sgal = load_phot(sgal, bands[i],  filter1=sfilters[i][0],
                                 filter2=sfilters[i][1], color_min=color_mins[i])

                # bin the data
                hist_attr = '%s_hist' % bands[i]
                bins_attr = '%s_bins' % bands[i]
                err_attr = '%s_err' % bands[i]
                mag_attr = '%s_mag' % bands[i]

                if bands[i] == 'ir':
                    # bin IRAC1, not IRAC2 (convention?)
                    mag_attr = '%s_mag1' % bands[i]


                for g in [gal, sgal]:
                    mag = g.__getattribute__(mag_attr)
                    hist, bins = np.histogram(mag, bins=bins)
                    err = np.sqrt(hist)

                    g.__setattr__(hist_attr, hist)
                    g.__setattr__(bins_attr, bins)
                    g.__setattr__(err_attr, err)

            # number of rgb and agb stars in the data
            for g in [gal, sgal]:
                g.nir_rgb, g.nir_agb = \
                    rgb_agb_regions(nir_offset, nir_trgb_exclude, nir_trgb,
                                    g.nir_mag)

                g.opt_rgb, g.opt_agb = \
                    rgb_agb_regions(opt_offset, opt_trgb_exclude, opt_trgb,
                                    g.opt_mag)

            # scale the simulation
            opt_norm, sc_opt_rgb, sc_opt_agb = \
                    normalize_simulation(sgal.opt_mag, float(len(gal.opt_rgb)),
                                         sgal.opt_rgb, sgal.opt_agb)
            nir_norm, sc_nir_rgb, sc_nir_agb = \
                    normalize_simulation(sgal.nir_mag, float(len(gal.nir_rgb)),
                                         sgal.nir_rgb, sgal.nir_agb)
            sgal.opt_hist *= opt_norm
            sgal.nir_hist *= nir_norm
            sgal.ir_hist *= nir_norm  # using the NIR scaling for IR.


            for i, band in enumerate(bands):
                j = 1  # choose filter2
                if band == 'ir':
                    j = 0  # choose IRAC1
                #hist_attr = '%s_hist' % band
                #bins_attr = '%s_bins' % band
                #err_attr = '%s_err' % band
                #ax = make_lf(gal.__getattribute__(hist_attr),
                #             gal.__getattribute__(bins_attr),
                #             sgal.__getattribute__(hist_attr),
                #             sgal.__getattribute__(bins_attr),
                #             gerrs=gal.__getattribute__(err_attr),
                #             serrs=sgal.__getattribute__(err_attr),
                #             ax=axs[i])
                ghist = np.histogram(gal.data[gfilters[i][j]][gal.icut], bins=bins)[0]
                shist = np.histogram(sgal.data[sfilters[i][j]][sgal.icut], bins=bins)[0]
                if band == 'opt':
                    shist *= opt_norm
                if 'ir' in band:
                    shist *= nir_norm
                gerr = np.sqrt(ghist)
                serr = np.sqrt(shist)
                ax = make_lf(ghist, bins, shist, bins, gerrs=gerr, serrs=serr,
                             ax=axs[i], colors=['k', colors[k]], label=agb_mod.replace('_pagb', ''))
                ax.set_xlabel(r'$%s$' % sfilters[i][j].replace('_', ''), fontsize=20)

                ax.set_xlim(24, 14)
                ax.set_ylim(1, 500)

                if band == 'nir':
                    ax = add_vlines_lf(ax, nir_offset, nir_trgb_exclude,
                                       nir_trgb)
                    #ax.set_xlim(20, 16)
                if band == 'opt':
                    ax = add_vlines_lf(ax, opt_offset, opt_trgb_exclude,
                                       opt_trgb)
                    ax.set_title(r'$\rm{B21F%02i}$' % field, fontsize=20)
                    #ax.set_xlim(24, 16)
                ax.legend(loc=0)
        plt.savefig('M31_B21_F%02i_lfs%s.png' % (field, extra))
        #[ax.set_ylim(1,300) for ax in [ax1, ax2]]
        '''
            coldats = ['CO', 'm_ini', 'MH']
            collabels = ['C/O', 'Initial Mass', '[M/H]']
            for coldat, collabel in zip(coldats, collabels):
                vmin = np.min([sgal.data[coldat][sgal.opt_agb].min(),
                               sgal.data[coldat][sgal.nir_agb].min()])
                vmax = np.max([sgal.data[coldat][sgal.opt_agb].max(),
                               sgal.data[coldat][sgal.nir_agb].max()])

                fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(12, 6))
                ax1 = gal.plot_cmd(gal.data[gopt_filt1][icut] -
                                    gal.data[gopt_filt2][icut],
                                   gal.data[gopt_filt2][icut], ax=ax1)
                ax2 = gal.plot_cmd(gal.data[gnir_filt1][icut] -
                                      gal.data[gnir_filt2][icut],
                                   gal.data[gnir_filt2][icut], ax=ax2)
                cs = ax1.scatter(sgal.opt_mag1[sgal.opt_agb] -
                                 sgal.opt_mag[sgal.opt_agb],
                                 sgal.opt_mag[sgal.opt_agb],
                                 c=sgal.data[coldat][sgal.opt_agb],
                                 cmap=plt.cm.get_cmap('RdBu'), vmin=vmin,
                                 vmax=vmax, alpha=0.3)

                cs = ax2.scatter(sgal.nir_mag1[sgal.nir_agb] -
                                    sgal.nir_mag[sgal.nir_agb],
                                 sgal.nir_mag[sgal.nir_agb],
                                 c=sgal.data[coldat][sgal.nir_agb],
                                 cmap=plt.cm.get_cmap('RdBu'),vmin=vmin,
                                 vmax=vmax,
                              alpha=0.3)
                cbar = plt.colorbar(cs)
                cbar.set_label(r'$%s$' % collabel)
                ax1.set_xlim(3.2, 5.5)
                ax1.set_ylim(20.5, 19.2)
                ax2.set_xlim(.95, 1.6)
                ax2.set_ylim(18.2, 15.1)
                ax1.set_xlabel(r'$%s-%s$' % (sopt_filt1, sopt_filt2),
                               fontsize=20)
                ax2.set_xlabel(r'$%s-%s$' % (snir_filt1, snir_filt2),
                               fontsize=20)
                ax1.set_ylabel(r'$%s$' % sopt_filt2, fontsize=20)
                ax2.set_ylabel(r'$%s$' % snir_filt2, fontsize=20)
                ax1.set_title('B21F%02i' % field)
                plt.savefig('M31_B21_F%02i_cmd_%s_%s%s.png' % (field, coldat,
                                                               agb_mod, extra))
        '''
# trgbs:
#ax1.vlines(nir_trgb, *ax1.get_ylim(), lw=2, color='gray')
#rgb_hist, rgb_bins = np.histogram(sgal.data.IRAC_36[sgal.data.stage==3.], bins=bins)
#ax2.plot(rgb_bins[1:], rgbhist, drawstyle=steps-mid, color='grey', lw=2)


# for using irac-phat:
# select data stars to compare with model
#ibest = list(set(gal.gst) & set(gal.nir_color_cut))# & set(gal.finds))

#gal.opt_mag = gal.data.F814W_VEGA[ibest]
#gal.opt_mag1 = gal.data.F475W_VEGA[ibest]
#gal.nir_mag1 = gal.data.F110W_VEGA[ibest]
#gal.nir_mag = gal.data.F160W_VEGA[ibest]
#gal.ir_mag1 = gal.data.IRAC1[ibest]
#gal.ir_mag2 = gal.data.IRAC2[ibest]


def plot_gal_3ways():
    ''' just put here so I don't have to delete it. Photometric issues tricked me'''
    bins = np.arange(10,28,0.1)
    gal.nir_hist, gal.nir_bins = np.histogram(gal.ir_mag1[gal.opt_agb], bins=bins)
    gal.ir_hist, gal.ir_bins = np.histogram(gal.ir_mag2[gal.opt_agb], bins=bins)
    nir_err = np.sqrt(gal.nir_hist)
    ir_err = np.sqrt(gal.ir_hist)
    ax1, ax2 = make_lf(gal.nir_hist, gal.nir_bins, gal.ir_hist, gal.ir_bins,
                       err1s=nir_err, err2s=ir_err, color='darkred',
                       label='opt cut')

    gal.nir_hist, gal.nir_bins = np.histogram(gal.ir_mag1[gal.nir_agb], bins=bins)
    gal.ir_hist, gal.ir_bins = np.histogram(gal.ir_mag2[gal.nir_agb], bins=bins)
    nir_err = np.sqrt(gal.nir_hist)
    ir_err = np.sqrt(gal.ir_hist)
    ax1, ax2 = make_lf(gal.nir_hist, gal.nir_bins, gal.ir_hist, gal.ir_bins,
                       err1s=nir_err, err2s=ir_err, color='navy', ax1=ax1, ax2=ax2,
                       label='nir cut')

    gal.nir_hist, gal.nir_bins = np.histogram(gal.ir_mag1, bins=bins)
    gal.ir_hist, gal.ir_bins = np.histogram(gal.ir_mag2, bins=bins)
    nir_err = np.sqrt(gal.nir_hist)
    ir_err = np.sqrt(gal.ir_hist)
    ax1, ax2 = make_lf(gal.nir_hist, gal.nir_bins, gal.ir_hist, gal.ir_bins,
                       err1s=nir_err, err2s=ir_err, color='black', ax1=ax1, ax2=ax2,
                       label='full')
    fig, axs = plt.subplots(ncols=2)
    inds1, = np.nonzero((gal.ir_mag1 < 17.) & (gal.ir_mag1 > 16.5))
    inds2, = np.nonzero((gal.ir_mag1 < 16.5) & (gal.ir_mag1 > 16.))

    inds3 = [i for i in inds1 if i in gal.nir_agb]
    inds4 = [i for i in inds2 if i in gal.nir_agb]

    for ax in axs:
        ax = gal.plot_cmd(gal.data.F110W_VEGA[gal.gst]-gal.data.F160W_VEGA[gal.gst],
                          gal.data.F160W_VEGA[gal.gst], ax=ax)

    axs[0].plot(gal.nir_mag1[inds4]-gal.nir_mag[inds4], gal.nir_mag[inds4], '.',
                color='orange', zorder=1000, alpha=0.5)

    axs[1].plot(gal.nir_mag1[inds3]-gal.nir_mag[inds3], gal.nir_mag[inds3], '.', color='navy', zorder=1000,
             alpha=0.5)
    for ax in axs:
        ax.set_xlim(0.5, 2)
        ax.set_ylim(21.5, 15)
        ax.set_xlabel(r'$F110W-F160W$', fontsize=20)
        ax.set_ylabel(r'$F160W$', fontsize=20)
