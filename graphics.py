import fileIO
import numpy as np
import os
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib import rcParams
from matplotlib.ticker import NullFormatter, MultipleLocator, ScalarFormatter
rcParams['text.usetex'] = True
rcParams['text.latex.unicode'] = False
rcParams['axes.linewidth'] = 1
rcParams['ytick.labelsize'] = 'large'
rcParams['xtick.labelsize'] = 'large'
rcParams['axes.edgecolor'] = 'black'
#rc('text', usetex=True)

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


def diag_plots(track, logl, logt, age, slopes, Qs, addpt, rising, fits,
               **kwargs):
    agb_mix = kwargs.get('agb_mix')
    set_name = kwargs.get('set_name')
    ext = '.png'
    logt_lim = (3.75, 3.35)
    logl_lim = (2.4, 4.8)
    lage_lim = (1., 1e7)
    co_lim = (0, 5)

    # HRD
    fig = plt.figure()
    ax = plt.axes()
    plotpath = os.path.join(kwargs.get('diagnostic_dir'), 'HRD/')
    fileIO.ensure_dir(plotpath)
    ax.annotate(r'$%s$' % agb_mix.replace('_', '\ '), xy=(3.43, 2.8))
    ax.annotate(r'$%s$' % set_name, xy=(3.43, 2.7))
    ax.annotate(r'$M=%.2f$' % track.mass, xy=(3.43, 2.6))
    ax.plot(logt, logl, color='black')

    ax.plot(logt[map(int, Qs)], logl[map(int, Qs)], color='green', lw=2)
    [ax.plot(logt[q], logl[q], 'o', color='green') for q in Qs]
    [ax.plot(logt[add], logl[add], 'o', color='purple')
     for add in addpt if len(addpt) > 0]
    ax.set_xlim(logt_lim)
    ax.set_ylim(logl_lim)
    ax.set_xlabel(r'$\log\ Te$')
    ax.set_ylabel(r'$\log\ L_{\odot}$')
    fname = os.path.split(track.name)[1].replace('.dat', '')
    fig_name = os.path.join(plotpath, '_'.join(('diag', fname)))
    plt.savefig('%s%s' % (fig_name, ext))
    plt.close()

    ax1, ax2, ax3 = set_up_three_panel_plot()
    # AGE v TE
    ax1.plot(age, logt, color='black')
    ax1.plot(age[Qs], logt[Qs], 'o', color='green')
    ax1.yaxis.set_major_locator(MultipleLocator(.1))
    ax1.yaxis.set_minor_locator(MultipleLocator(.05))

    # AGE v L
    ax2.plot(age, logl, color='black')
    ax2.plot(age[Qs], logl[Qs], 'o', color='green')
    ax2.yaxis.set_major_locator(MultipleLocator(.2))
    ax2.yaxis.set_minor_locator(MultipleLocator(.1))

    # AGE vs CO
    CO = track.get_col('CO')
    ax3.plot(age, CO, color='black')
    ax3.plot(age[Qs], CO[Qs], 'o', color='green')
    majorFormatter = ScalarFormatter()
    majorFormatter.set_powerlimits((-3, 4))
    ax1.xaxis.set_major_formatter(majorFormatter)
    ax3.xaxis.set_major_formatter(nullfmt)
    ax2.xaxis.set_major_formatter(nullfmt)
    ax3.annotate(r'\rm{%s}' % agb_mix.replace('_', '\ '), xy=(.06, .87),
                 textcoords =('axes fraction'))
    ax3.annotate(r'\rm{%s}' % set_name.replace('_', '\ '), xy=(.06, .77),
                 textcoords =('axes fraction'))
    ax3.annotate(r'$M=%.2f$' % track.mass, xy=(.06, .67),
                 textcoords =('axes fraction'))
    ax3.set_ylabel(r'$C/O$')
    ax2.set_ylabel(r'$\log\ L_{\odot}$')
    ax1.set_xlabel(r'$\rm{Age (yr)}$')
    ax1.set_ylabel(r'$\log\ Te$')

    plotpath = os.path.join(kwargs.get('diagnostic_dir'), 'age_v/')
    fileIO.ensure_dir(plotpath)
    fname = os.path.split(track.name)[1].replace('.dat', '')
    fig_name = os.path.join(plotpath, '_'.join(('diag', fname)))
    
    majorLocator = MultipleLocator()
    for ax in ax1, ax2, ax3:
        ax.grid(linestyle=': ', color='grey')

    plt.savefig('%s%s%s' % (fig_name, '_age_v', ext))
    plt.close()
    return


def plot_ifmr(imfrfile):
    mi, mf, z = np.loadtxt(imfrfile, unpack=True)
    zinds, unz = get_unique_inds(z)
    cols = discrete_colors(len(zinds))
    [plt.plot(mi[zinds[i]], mf[zinds[i]], color=cols[i], label=str(z[unz[i]]))
        for i in range(len(unz))]
    [plt.plot(mi[zinds[i]], mf[zinds[i]], 'o', ms=4, color=cols[i],
              mec='white', alpha=0.5) for i in range(len(unz))]
    plt.legend(loc=2, frameon=False)
    plt.xlabel(r'$M_i/M_{\odot}$', fontsize=20)
    plt.ylabel(r'$M_f/M_{\odot}$', fontsize=20)
    plt.savefig(imfrfile.replace('dat', 'png'))
    plt.close()
    return


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

    
def plot_cluster_test(lifetimesfile, **kwargs):
    agbtrack_dir = kwargs.get('agbtrack_dir')
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
