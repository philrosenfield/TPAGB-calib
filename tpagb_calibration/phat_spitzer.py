import ResolvedStellarPops as rsp
import pyfits
import numpy as np
import matplotlib.pylab as plt
from pop_synth.stellar_pops import normalize_simulation, rgb_agb_regions
plt_labelsize = 20
tck_labelsize = 16


def load_data(fname, opt_color_min, nir_color_min, field=None, ir_filter1=None,
              ir_filter2=None):
    gal = rsp.StarPop()
    gal.data = pyfits.getdata(fname)

    gst, = np.nonzero((gal.data.F160W_GST == 'T') & (gal.data.F110W_GST == 'T') &
                      (gal.data.F475W_GST == 'T') & (gal.data.F814W_GST == 'T'))
    gal.gst = gst

    if field is not None:
        finds, = np.nonzero(gal.data.field('field') == field)
        # If in overlapping fields, I added the field flag together.
        if field == 6 or field == 12:
            finds2, = np.nonzero(gal.data.field('field') == 18)
            finds = np.concatenate([finds,finds2])
        gal.finds = finds

    gal = load_phot(gal, 'F475W_VEGA', 'F814W_VEGA', 'F110W_VEGA', 'F160W_VEGA',
                    opt_color_min, nir_color_min, ir_filter1=ir_filter1,
                    ir_filter2=ir_filter2)
    return gal


def load_phot(StarPop, opt_filter1, opt_filter2, nir_filter1, nir_filter2,
              opt_color_min, nir_color_min, ir_filter1=None, ir_filter2=None ):
    StarPop.opt_color = StarPop.data[opt_filter1] - StarPop.data[opt_filter2]
    StarPop.nir_color = StarPop.data[nir_filter1] - StarPop.data[nir_filter2]
    if not None in [ir_filter1, ir_filter2]:
        StarPop.ir_color = StarPop.data[ir_filter1] - StarPop.data[ir_filter2]

    StarPop.opt_color_cut, = np.nonzero((StarPop.opt_color) > opt_color_min)
    StarPop.nir_color_cut, = np.nonzero((StarPop.nir_color) > nir_color_min)
    return StarPop


def load_model(fname, opt_color_min, nir_color_min):
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
                              names=cols).view(np.recarray)
    sgal = load_phot(sgal, 'F475W', 'F814W', 'F110W', 'F160W', opt_color_min,
                     nir_color_min, ir_filter1='IRAC_36', ir_filter2='IRAC_45')
    return sgal


def make_lf(hist1, bin1s, hist2, bin2s, err1s=None, err2s=None, ax1=None,
            ax2=None, color='black', label=''):

    if ax1 is None and ax2 is None:
        fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(12,6))
    if not None in [err1s, err2s]:
        plt_kw = {'drawstyle': 'steps-mid', 'color': color,
                  'lw': 2, 'label': label}
        ax1.errorbar(bin1s[1:], hist1, yerr=err1s, **plt_kw)
        ax2.errorbar(bin2s[1:], hist2, yerr=err2s, **plt_kw)
    else:
        plt_kw = {'linestyle': 'steps-mid', 'color': color,
                  'lw': 2, 'label': label}
        ax1.plot(bin1s[1:], hist1, **plt_kw)
        ax2.plot(bin2s[1:], hist2, **plt_kw)

    for ax in [ax1, ax2]:
        ax.set_yscale('log')
        ax.tick_params(labelsize=16)

    ax1.set_ylabel(r'${\rm Number\ of\ Stars}$', fontsize=20)
    return ax1, ax2

def add_vlines_lf(ax1, ax2, offsets, trgb_excludes, opt_trgb, nir_trgb):
    trgb = [opt_trgb, nir_trgb]
    for i, ax  in enumerate([ax1, ax2]):
        yarr = np.linspace(*ax.get_ylim())
        ax.fill_betweenx(yarr, trgb[i] - trgb_excludes[i],
                         trgb[i] + trgb_excludes[i],
                         color='black', alpha=0.1)
        ax.vlines(trgb[i], *ax.get_ylim(), color='black', linestyle='--')
        ax.vlines(offsets[i], *ax.get_ylim(), color='black', linestyle='--')
    return ax1, ax2

offsets = [22.5, 19.2]
trgb_excludes = [0.1, 0.2]
#trgb_excludes = [0., 0.]
opt_trgb = 20.5
nir_trgb = 18.3
opt_color_min = 3.5
nir_color_min = 0.98

# load data and model
#gal_file = '/home/rosenfield/research/TP-AGBcalib/PHAT/data/Spitzer/B21_phat_irac_fields.fits'
#ir_filt1 = 'IRAC1'
#ir_filt2 = 'IRAC2'
gal_file = '/home/rosenfield/research/TP-AGBcalib/PHAT/data/Spitzer/b21-6filt-cut-shallow-fields.fits'
ir_filt1 = None
ir_filt2 = None

fields = [6, 12, 15]
for field in fields:
    gal = load_data(gal_file, opt_color_min, nir_color_min, ir_filter1=ir_filt1,
                    ir_filter2=ir_filt2, field=field)

    gopt_filt1 = 'F475W_VEGA'
    gopt_filt2 = 'F814W_VEGA'
    gnir_filt1 = 'F110W_VEGA'
    gnir_filt2 = 'F160W_VEGA'

    # select data stars to compare with model
    icut = list(set(gal.gst) & set(gal.finds))
    opt_cut = list(set(icut) & set(gal.opt_color_cut))
    nir_cut = list(set(icut) & set(gal.nir_color_cut))
    gal.opt_mag1 = gal.data[gopt_filt1][opt_cut]
    gal.opt_mag = gal.data[gopt_filt2][opt_cut]
    gal.nir_mag1 = gal.data[gnir_filt1][nir_cut]
    gal.nir_mag = gal.data[gnir_filt2][nir_cut]

    # number of rgb and agb stars in the data
    gal.opt_rgb, gal.nir_rgb, gal.opt_agb, gal.nir_agb = \
        rgb_agb_regions(offsets, trgb_excludes, opt_trgb, nir_trgb,
                        gal.opt_mag, gal.nir_mag)

    # bin the data
    bins = np.arange(16, 28, 0.1)
    gal.opt_hist, gal.opt_bins = np.histogram(gal.opt_mag, bins=bins)
    gal.nir_hist, gal.nir_bins = np.histogram(gal.nir_mag, bins=bins)
    opt_err = np.sqrt(gal.opt_hist)
    nir_err = np.sqrt(gal.nir_hist)

    # the best sfh with aringer
    #sgal_name = '/home/rosenfield/research/TP-AGBcalib/PHAT/vary_sfh/output_M31-B21_3x6-006.aringer.dat'
    # the best sfh
    sgal_name = '/home/rosenfield/research/TP-AGBcalib/PHAT/vary_sfh/output_M31-B21_3x6-0%02i.dat' % field

    sgal = load_model(sgal_name, opt_color_min, nir_color_min)
    sopt_filt1 = 'F475W'
    sopt_filt2 = 'F814W'
    snir_filt1 = 'F110W'
    snir_filt2 = 'F160W'

    sgal.opt_mag1 = sgal.data[sopt_filt1][sgal.opt_color_cut]
    sgal.opt_mag = sgal.data[sopt_filt2][sgal.opt_color_cut]
    sgal.nir_mag1 = sgal.data[snir_filt1][sgal.nir_color_cut]
    sgal.nir_mag = sgal.data[snir_filt2][sgal.nir_color_cut]

    sgal.opt_hist, sgal.opt_bins = np.histogram(sgal.opt_mag, bins=bins)
    sgal.nir_hist, sgal.nir_bins = np.histogram(sgal.nir_mag, bins=bins)

    sgal.opt_rgb, sgal.nir_rgb, sgal.opt_agb, sgal.nir_agb = \
        rgb_agb_regions(offsets, trgb_excludes, opt_trgb, nir_trgb,
                        sgal.opt_mag, sgal.nir_mag)

    # scale the simulation
    opt_norm, nir_norm, sc_opt_rgb, sc_nir_rgb, sc_opt_agb, sc_nir_agb = \
        normalize_simulation(sgal.opt_mag, sgal.nir_mag, float(len(gal.opt_rgb)),
                             float(len(gal.nir_rgb)), sgal.opt_rgb, sgal.nir_rgb,
                             sgal.opt_agb, sgal.nir_agb)

    sgal.opt_hist *= opt_norm
    sgal.nir_hist *= nir_norm

    ax1, ax2 = make_lf(gal.opt_hist, gal.opt_bins, gal.nir_hist, gal.nir_bins,
                           err1s=opt_err, err2s=nir_err)

    ax1, ax2 = make_lf(sgal.opt_hist, sgal.opt_bins, sgal.nir_hist, sgal.nir_bins,
                       color='darkred', label='sim', ax1=ax1, ax2=ax2)

    ax1.set_xlabel(r'$%s$' % sopt_filt2, fontsize=20)
    ax2.set_xlabel(r'$%s$' % snir_filt2, fontsize=20)
    ax1.set_xlim(23, 16)
    ax2.set_xlim(21.3, 14)
    [ax.set_ylim(1, ax.get_ylim()[1]) for ax in [ax1, ax2]]
    ax1, ax2 = add_vlines_lf(ax1, ax2, offsets, trgb_excludes, opt_trgb, nir_trgb)
    ax1.set_title('B21F%02i' % field)
    plt.savefig('M31_B21_F%02i_lfs.png' % field)
    #[ax.set_ylim(1,300) for ax in [ax1, ax2]]

    coldats = ['CO', 'm_ini', 'MH']
    collabels = ['C/O', 'Initial Mass', '[M/H]']
    for coldat, collabel in zip(coldats, collabels):
        vmin = np.min([sgal.data[coldat][sgal.opt_agb].min(),
                       sgal.data[coldat][sgal.nir_agb].min()])
        vmax = np.max([sgal.data[coldat][sgal.opt_agb].max(),
                       sgal.data[coldat][sgal.nir_agb].max()])

        fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(12, 6))
        ax1 = gal.plot_cmd(gal.data[gopt_filt1][icut] - gal.data[gopt_filt2][icut],
                           gal.data[gopt_filt2][icut], ax=ax1)
        ax2 = gal.plot_cmd(gal.data[gnir_filt1][icut] - gal.data[gnir_filt2][icut],
                           gal.data[gnir_filt2][icut], ax=ax2)
        cs = ax1.scatter(sgal.opt_mag1[sgal.opt_agb] - sgal.opt_mag[sgal.opt_agb],
                         sgal.opt_mag[sgal.opt_agb],
                         c=sgal.data[coldat][sgal.opt_agb],
                         cmap=plt.cm.get_cmap('RdBu'), vmin=vmin, vmax=vmax,
                         alpha=0.3)

        cs = ax2.scatter(sgal.nir_mag1[sgal.nir_agb] - sgal.nir_mag[sgal.nir_agb],
                         sgal.nir_mag[sgal.nir_agb],
                         c=sgal.data[coldat][sgal.nir_agb],
                         cmap=plt.cm.get_cmap('RdBu'),vmin=vmin, vmax=vmax,
                      alpha=0.3)
        cbar = plt.colorbar(cs)
        cbar.set_label(r'$%s$' % collabel)
        ax1.set_xlim(3.2, 5.5)
        ax1.set_ylim(20.5, 19.2)
        ax2.set_xlim(.95, 1.6)
        ax2.set_ylim(18.2, 15.1)
        ax1.set_xlabel(r'$%s-%s$' % (sopt_filt1, sopt_filt2), fontsize=20)
        ax2.set_xlabel(r'$%s-%s$' % (snir_filt1, snir_filt2), fontsize=20)
        ax1.set_ylabel(r'$%s$' % sopt_filt2, fontsize=20)
        ax2.set_ylabel(r'$%s$' % snir_filt2, fontsize=20)
        ax1.set_title('B21F%02i' % field)
        plt.savefig('M31_B21_F%02i_cmd_%s.png' % (field, coldat))

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
