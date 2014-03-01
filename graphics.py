import fileIO
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rcParams
from matplotlib.ticker import NullFormatter, MultipleLocator, ScalarFormatter, MaxNLocator
rcParams['text.usetex'] = True
rcParams['text.latex.unicode'] = False
rcParams['axes.linewidth'] = 1
rcParams['ytick.labelsize'] = 'large'
rcParams['xtick.labelsize'] = 'large'
rcParams['axes.edgecolor'] = 'black'
#rc('text', usetex=True)
import ResolvedStellarPops.graphics.GraphicsUtils as rspg

nullfmt = NullFormatter()  # no labels

def discrete_colors(ncolors, colormap='gist_rainbow'):
    colors = []
    cmap = cm.get_cmap(colormap)
    # color will now be an RGBA tuple
    colors = [(cmap(1. * i / ncolors)) for i in range(ncolors)]
    return colors


def set_up_three_panel_plot():
    fig = plt.figure(figsize=(8, 8))

    d = 0.02
    low_b = 0.1

    left, width = 0.1, 0.8
    dh = (1. - 3 * d - low_b - .05) / 3.

    low_t = low_b + dh
    mid_b = low_t + d
    mid_t = mid_b + dh
    top_b = mid_t + d
    top_t = top_b + dh

    bottoms = [low_b, mid_b, top_b]

    axs = [plt.axes([left, b, width, dh]) for b in bottoms]
    return axs


def two_panel_plot_vert():
    fig = plt.figure(2, figsize=(8, 8))
    left, width = 0.13, 0.83
    bottom, height = 0.1, 0.41
    dh = 0.03

    axis1 = [left, bottom, width, height]
    axis2 = [left, (bottom + height + dh), width, height]

    ax1 = plt.axes(axis1)
    ax2 = plt.axes(axis2)
    ax2.xaxis.set_major_formatter(nullfmt)

    return ax1, ax2


def diag_plots(track, infile):
    agb_mix = infile.agb_mix
    set_name = infile.set_name
    ext = '.png'
    logt_lim = (3.75, 3.35)
    logl_lim = (2.4, 4.8)
    lage_lim = (1., 1e7)
    co_lim = (0, 5)

    logl = track.get_col('L_star')
    logt = track.get_col('T_star')
    addpt = track.addpt
    Qs = list(track.Qs)
    # HRD
    fig = plt.figure()
    ax = plt.axes()
    plotpath = os.path.join(infile.diagnostic_dir, 'HRD/')
    fileIO.ensure_dir(plotpath)
    ax.annotate(r'$%s$' % agb_mix.replace('_', '\ '), xy=(3.43, 2.8))
    ax.annotate(r'$%s$' % set_name, xy=(3.43, 2.7))
    ax.annotate(r'$M=%.2f$' % track.mass, xy=(3.43, 2.6))
    ax.plot(logt, logl, color='black')

    ax.plot(logt[Qs], logl[Qs], color='green', lw=2)
    ax.plot(logt[Qs], logl[Qs], 'o', color='green')
    if len(addpt) > 0:
        ax.plot(logt[addpt], logl[addpt], 'o', color='purple')
    ax.set_xlim(logt_lim)
    ax.set_ylim(logl_lim)
    ax.set_xlabel(r'$\log\ Te$')
    ax.set_ylabel(r'$\log\ L_{\odot}$')
    fname = os.path.split(track.name)[1].replace('.dat', '')
    fig_name = os.path.join(plotpath, '_'.join(('diag', fname)))
    plt.savefig('%s%s' % (fig_name, ext))
    plt.close()

    fig, (axs) = plt.subplots(nrows=3)
    ycols = ['logl', 'logt', 'C/O']
    annotate = [False, False, True]
    save_plot = [False, False, True]
    xlabels = [False, False, True]
    for i in range(len(ycols)):
        age_vs_plot(track, infile, ycol=ycols[i], ax=axs[i], annotate=annotate[i],
                    xlabels=xlabels[i], ylabels=True, save_plot=save_plot[i])


def age_vs_plot(track, infile, ycol='logl', ax=None, annotate=True, xlabels=True,
                save_plot=True, ylabels=True):
    agb_mix = infile.agb_mix
    set_name = infile.set_name

    if ycol == 'logl':
        ydata= track.get_col('L_star')
        majL = MultipleLocator(.2)
        minL = MultipleLocator(.1)
        ylab = '$\log\ L_{\odot}$'
    elif ycol == 'logt':
        ydata = track.get_col('T_star')
        majL = MultipleLocator(.1)
        minL = MultipleLocator(.05)
        ylab = '$\log\ Te$'
    elif ycol == 'C/O':
        ydata = track.get_col('CO')
        majL = MaxNLocator(4)
        minL = MaxNLocator(2)
        ylab = '$C/O$'
    else:
        print 'logl, logt, C/O only choices for y.'
        return

    age = track.get_col('ageyr')
    addpt = track.addpt
    Qs = list(track.Qs)

    if ax is None:
        fig, ax = plt.subplots()
    ax.plot(age, ydata, color='black')
    ax.plot(age[Qs], ydata[Qs], 'o', color='green')
    if len(addpt) > 0:
        ax.plot(age[addpt], ydata[addpt], 'o', color='purple')
    ax.yaxis.set_major_locator(majL)
    ax.yaxis.set_minor_locator(minL)
    majorFormatter = ScalarFormatter()
    majorFormatter.set_powerlimits((-3, 4))
    ax.xaxis.set_major_formatter(majorFormatter)

    if annotate is True:
        ax.text(0.06, 0.87, '${\\rm %s}$' % agb_mix.replace('_', '\ '),
                transform=ax.transAxes)
        ax.text(0.06, 0.77,'${\\rm %s}$' % set_name.replace('_', '\ '),
                transform=ax.transAxes)
        ax.text(0.06, 0.67, '$M=%.2f$' % track.mass,
                transform=ax.transAxes)
    if ylabels is True:
        ax.set_ylabel('$%s$' % ylab, fontsize=20)
    if xlabels is True:
        ax.set_xlabel('$\rm{Age (yr)}$', fontsize=20)

    if save_plot is True:
        plotpath = os.path.join(infile.diagnostic_dir, 'age_v/')
        fileIO.ensure_dir(plotpath)
        fname = os.path.split(track.name)[1].replace('.dat', '')
        fig_name = os.path.join(plotpath, '_'.join(('diag', fname)))
        plt.savefig('%s_age_v.png' % fig_name, dpi=300)
        plt.close()
    return

def bigplots(agb_tracks, infile):
    if type(agb_tracks[0]) == str:
        agb_tracks = [fileIO.get_numeric_data(a) for a in agb_tracks]
    out_fig = os.path.join(infile.diagnostic_dir, 'hrd_%.4f.png' % agb_tracks[0].metallicity)
    plot_title = '$\dot{M}_{\\rm M13}\ Z=%.4f$' % agb_tracks[0].metallicity
    nagb_tracks = len(agb_tracks)
    fig, (axs) = rspg.setup_multiplot(nagb_tracks,
                                      ylabel='$\log\ L\ (L_\odot)$',
                                      xlabel='$\log\ T_{\\rm eff}\ (K)$',
                                      title =plot_title,
                                      subplots_kwargs={'figsize': (30,30),
                                                       'squeeze': True})
    axs = axs.flatten()
    [hrd_slopes(agb_tracks[i], ax=axs[i]) for i in range(nagb_tracks)]
    plt.savefig(out_fig, dpi=300)
    ylabel = ['$C/O$', '$\log\ L\ (L_\odot)$', '$\log\ T_{\\rm eff}\ (K)$']
    ycol = ['C/O', 'logl', 'logt']
    for i in range(len(ylabel)):
        fig, (axs) = rspg.setup_multiplot(nagb_tracks,
                                          ylabel=ylabel[i],
                                          xlabel='${\\rm Age (yr)}$',
                                          title =plot_title,
                                          subplots_kwargs={'figsize': (30,30),
                                                           'squeeze': True})
        axs = axs.flatten()
        [age_vs_plot(agb_tracks[j], infile, ycol=ycol[i], ax=axs[j], annotate=True,
                     xlabels=False, ylabels=False, save_plot=False)
         for j in range(len(agb_tracks))]
        out_fig = out_fig.replace('hrd', 'age_v_%s' % ycol[i].lower().replace('c/o',
                                                                              'co'))
        plt.savefig(out_fig, dpi=300)

def hrd_slopes(track, ax=None):
    logl = track.get_col('L_star')
    logt = track.get_col('T_star')
    Qs = list(track.Qs)
    x = np.linspace(min(logt), max(logt))
    if ax is None:
        fig, ax = plt.subplots()
    [ax.plot(x, track.fits[i][0]*x + track.fits[i][1], '--', color='blue',
             lw=0.5) for i in range(len(track.fits))]

    ax.plot(logt, logl, color='black')
    [ax.plot(np.sort(logt[r]), logl[r][np.argsort(logt[r])], color='red', lw=2) for r in track.rising
     if len(r) > 0]
    ax.plot(logt[Qs], logl[Qs], color='green', lw=2)
    [ax.plot(logt[q], logl[q], 'o', color='green') for q in Qs]
    [ax.plot(logt[add], logl[add], 'o', color='purple')
     for add in track.addpt if len(track.addpt) > 0]
    ax.axis([max(logt) + 0.01, min(logt) - 0.01, min(logl) - 0.01,
             max(logl) + 0.01])
    ax.text(0.1, 0.85, '$%.3f$' % track.mass, transform=ax.transAxes, fontsize=16)

def check_output():
    pass

def plot_ifmr(imfrfile, ax=None, zs=None, data_mi=None, data_mf=None, data_mierr=None, data_mferr=None):
    mi, mf, z = np.loadtxt(imfrfile, unpack=True)
    zinds, unz = get_unique_inds(z)
    if zs is None:
        zs = np.unique(z)
    cols = discrete_colors(len(zs))
    cols = cols * len(zinds)
    if ax is None:
        fig, ax = plt.subplots()
    [ax.plot(mi[zinds[i]], mf[zinds[i]], color=cols[i], label=str(z[unz[i]]))
        for i in range(len(unz)) if z[unz[i]] in zs]
    [ax.plot(mi[zinds[i]], mf[zinds[i]], 'o', ms=4, color=cols[i],
              mec='white', alpha=0.5) for i in range(len(unz)) if z[unz[i]] in zs]
    if data_mf is not None:
        if data_mierr is not None:
            from matplotlib.patches import Rectangle
            bottom_left = (data_mi - data_mierr/2., data_mf - data_mferr/2.)
            width = data_mierr
            height = data_mferr
            rect = Rectangle(bottom_left, width, height, color='black', alpha=.25)
            ax.add_patch(rect)
    ax.legend(loc=2, frameon=False)
    ax.set_xlabel(r'$M_i/M_{\odot}$', fontsize=20)
    ax.set_ylabel(r'$M_f/M_{\odot}$', fontsize=20)
    plt.savefig(imfrfile.replace('dat', 'png'))
    #plt.close()
    return ax


def two_panel_plot_vert(fign=2):
    fig = plt.figure(fign, figsize=(8, 8))
    left, width = 0.13, 0.83
    bottom, height = 0.1, 0.41
    dh = 0.03

    axis1 = [left, bottom, width, height]
    axis2 = [left, (bottom + height + dh), width, height]

    ax1 = plt.axes(axis1)
    ax2 = plt.axes(axis2)
    ax2.xaxis.set_major_formatter(nullfmt)
    return ax1, ax2


def plot_cluster_test(lifetimesfile, infile):
    agbtrack_dir = infile.agbtrack_dir
    cluster_data = os.path.join(agbtrack_dir, 'cluster_data.dat')
    cmass, ctauc, ecp, ecm, ctaum, emp, emm = np.loadtxt(cluster_data,
                                                         unpack=True)
    z, mass, tauc, taum = np.loadtxt(lifetimesfile, unpack=True)
    zinds, unz = get_unique_inds(z)
    cols = discrete_colors(len(zinds))

    ax1, ax2 = two_panel_plot_vert(fign=3)
    ax2.set_ylabel(r'$\tau_M\ (\rm{Myr})$', fontsize=20)
    ax1.set_ylabel(r'$\tau_C\ (\rm{Myr})$', fontsize=20)
    ax1.set_xlabel(r'$M_{TO}\ (\rm{M}_\odot)$', fontsize=20)
    ax1.errorbar(cmass, ctauc, yerr=[ecm, ecp], fmt='o', color='black',
                 label='LMC')
    ax2.errorbar(cmass, ctaum, yerr=[emm, emp], fmt='o', color='black')

    for i in range(len(unz)):
        ax1.plot(mass[zinds[i]], tauc[zinds[i]], color=cols[i],
                 label=str(z[unz[i]]))
        ax2.plot(mass[zinds[i]], taum[zinds[i]], color=cols[i])

    ax1.legend(loc=0, frameon=False, numpoints=1)
    plt.savefig(lifetimesfile.replace('dat', 'png'))
    print 'wrote', lifetimesfile.replace('dat', 'png')
    plt.close()
    return

def get_unique_inds(ntp):
    un = np.unique(ntp)
    if un.size == 1:
        print 'only one themal pulse.'
        return un, list(ntp).index(u)
    else:
        # this is the first step in each TP.
        iTPs = [list(ntp).index(u) for u in un]
        # The indices of each TP.
        TPs = [np.arange(iTPs[i], iTPs[i + 1]) for i in range(len(iTPs) - 1)]
        # don't forget the last one.
        TPs.append(np.arange(iTPs[i + 1], len(ntp)))
        return TPs, iTPs
