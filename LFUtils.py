#
#  LumnFncfromCMD.py
#
#  Created by Philip Rosenfield on 11/22/10.
#
import os
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.patheffects import withStroke
from matplotlib import rc
myeffect = withStroke(foreground="w", linewidth=3)
kwargs = dict(path_effects=[myeffect])
rc('text', usetex=True)
#from ResolvedStellarPops.math_utils import bayesian_blocks

# Phil's functions:
from TPAGBparams import *


def read_mettable():
    tab = os.path.join(table_src, 'IR_NAGBs.dat')
    dtype = [('fitstable', '|S46'),
             ('IR_TRGB', '<f8'),


             ('N_AGB', '<f8'),
             ('logOH', '<f8'),
             ('OHerr', '<f8'),
             ('Z', '<f8'),
             ('FeH', '<f8'),
             ('Ttype', '<f8')]
    table = np.genfromtxt(tab, dtype=dtype, delimiter=',')
    return table


def get_key_fromtable(ID, key):
    tab = read_mettable()
    names = list(tab['fitstable'])
    i, = [names.index(i) for i in names if ID in i]
    return tab[key][i]


def get_trgb_ir_nAGB(Target):
    trgb_ir = get_key_fromtable(Target, 'IR_TRGB')
    nAGB = get_key_fromtable(Target, 'N_AGB')
    return trgb_ir, nAGB


def get_chi2(obs, exp):
    return np.sum((obs - exp) ** 2 / exp)


def setup_lfplot(gal):
    # plot limits
    bottom = 0.1
    height = 0.8
    widths = [0.28, 0.28, 0.2]
    lefts = [0.1, 0.41, 0.72]

    fig = plt.figure(1, figsize=(9, 9))

    axs = [plt.axes([lefts[i], bottom, widths[i], height])
           for i in range(3)]
    lab_kw = {'fontsize': 20}
    # titles
    axs[0].set_title(r'$Data$', color='black', **lab_kw)
    axs[1].set_title(r'$Model$', color='red', **lab_kw)
    axs[0].set_xlabel(r'$%s-%s$' % (gal.filter1, gal.filter2), **lab_kw)
    axs[0].set_ylabel(r'$%s$' % gal.filter2, **lab_kw)
    axs[1].set_xlabel(axs[0].get_xlabel(), **lab_kw)
    axs[2].set_xlabel(r'$\#$', **lab_kw)
    # no formatters on mid and right plots
    [ax.yaxis.set_major_formatter(NullFormatter()) for ax in axs[1:]]

    for ax in axs:
        ax.tick_params(labelsize=16)
    return (fig, axs)


def plot_lines(axs, xrange, yval):
    [ax.plot((xrange), (yval, yval), color='black') for ax in axs]
    return


def plot_numbs(ax, item, xpos, ypos, **kwargs):
    ax.annotate(r'$%i$' % item, xy=(xpos, ypos), ha='left',
                **kwargs)
    return


def make_title(gal, fig):
    text_kwargs = {'ha': 'center', 'va': 'top', 'size': 20}
    if np.isfinite(gal.z):
        title = r'$\rm{%s\ m}_{TRGB}=%.3f\ Z=%.4f$' % (gal.target, gal.trgb,
                                                       gal.z)
    else:
        title = r'$\rm{%s\ m}_{TRGB}=%.3f\ Z=...$' % (gal.target, gal.trgb)

    fig.text(0.5, 0.96, title, **text_kwargs)
    return


def setup_plot_numbers():
    dim, bright = maglims

    sgal.load_ic_mstar()

    bright_limit = gal.trgb - sgal.count_offset

    bright_mstars = rsp.math_utils.brighter(sgal.ast_mag2, bright_limit,
                                            inds=sgal.imstars)

    bright_cstars = rsp.math_utils.brighter(sgal.ast_mag2, bright_limit,
                                            inds=sgal.imcstars)

    sim_stars = sgal.ast_mag2[sgal.rec][sgal.rgb_norm_inds]
    nbright_rgb = sgal.rel_agb.size - len(bright_mstars) - len(bright_cstars)

    # model used for normalization
    norm_inds = list(set(sgal.norm) & set(sgal.rel_ind))

    if sgal.rel_agb.size - len(bright_mstars) - len(bright_cstars) - nbright_rgb != 0.:
        print ''
        print sgal.rel_agb.size, len(bright_mstars) - len(bright_cstars) - nbright_rgb
        print ''


def annotate_line(ax, yval, str, offset=0.1, text_kw={}):
    text_kw = dict({'fontsize': 20}.items() + text_kw.items())
    ax.text(ax.get_xlim()[0]+offset, yval-offset, str, **text_kw)


def plot_LFIR(gal, sgal, p_value, maglims, res=0.1):
    model_color = 'red'
    data_color = 'black'
    nbins = (gal.mag2.max() - gal.mag2.min()) / res
    gal_hist, bins = np.histogram(gal.mag2, nbins)
    sgal_hist, _ = np.histogram(sgal.ast_mag2[sgal.rec], bins=bins)
    sgal_hist *= sgal.rgb_norm
    mhist, _ = np.histogram(sgal.ast_mag2[sgal.rec][sgal.imstar], nbins)
    chist, _ = np.histogram(sgal.ast_mag2[sgal.rec][sgal.icstar], bins=bins)

    (fig, [axs]) = setup_lfplot(gal)

    make_title(gal, fig)
    plt_kw = {'threshold': 25, 'levels': 3, 'scatter_args': {'alpha': 1}}
    # plot data
    gal.plot_cmd(gal.color, gal.mag2, ax=axs[0], color=data_color, **plt_kw)
    # plot simulation
    sgal.plot_cmd(sgal.ast_color[sgal.rec][sgal.rgb_norm_inds],
                  sgal.ast_mag2[sgal.rec][sgal.rgb_norm_inds], ax=axs[1],
                  color=model_color, **plt_kw)
    # plot histogram
    axs[2].semilogx(gal_hist, bins[1:], drawstyle='steps', color=data_color, lw=2)
    axs[2].semilogx(sgal_hist, bins[1:], drawstyle='steps', color=model_color, lw=2)

    # fix axes
    axs[0].set_ylim(axs[0].get_ylim()[0], 20)
    [ax.set_ylim(axs[0].get_ylim()) for ax in axs[1:]]
    axs[1].set_xlim(axs[0].get_xlim())

    # lines and numbers on plots
    data_agb = rsp.math_utils.brighter(gal.mag2, gal.trgb).size
    line_on_it_kw = {'annotate': 0, 'ls': '-'}
    gal.put_a_line_on_it(axs[0], gal.trgb, color=data_color, **line_on_it_kw)
    annotate_line(axs[0], gal.trgb, '%i' % data_agb, text_kw={'color': data_color})

    sim_agb = rsp.math_utils.brighter(sgal.ast_mag2[sgal.rec][sgal.rgb_norm_inds], gal.trgb).size
    gal.put_a_line_on_it(axs[1], gal.trgb, color=model_color, **line_on_it_kw)
    annotate_line(axs[1], gal.trgb, '%i' % sim_agb, text_kw={'color': model_color})
    
    '''
    axm.plot(sgal.ast_color[mstars], sgal.ast_mag2[mstars], '.',
             color='darkblue', mew=0, label=r'$M=%i$' % (len(bright_mstars)))

    axm.plot(sgal.ast_color[cstars], sgal.ast_mag2[cstars], '.',
             color='darkblue', mew=0, label=r'$C=%i$' % (len(bright_cstars)))

    axm.legend(frameon=False, loc=2, numpoints=1)

    #[plot_lines([axd, axm], axd.get_xlim(), m) for m in (dim, bright)]
    # also?
    plot_lines([axd, axm], axd.get_xlim(), gal.trgb)
    plot_lines([axd, axm], axd.get_xlim(), gal.trgb-sgal.count_offset)

    kwargs['color'] = 'black'
    text_offset = 0.02
    xpos = axd.get_xlim()[0] + 2 * text_offset
    yposs = np.asarray(maglims) - text_offset
    ypos = gal.trgb-text_offset
    plot_numbs(axd, gal.iagb.size-sgal.count_offset, xpos, ypos, **kwargs)
    plot_numbs(axd, gal.rgb_norm.size, xpos, yposs[0], **kwargs)

    kwargs['color'] = 'red'
    xpos = axm.get_xlim()[0] + text_offset
    plot_numbs(axm, sgal.rel_agb.size, xpos, ypos, **kwargs)
    plot_numbs(axm, len(norm_inds), xpos, yposs[0], **kwargs)

    # plot hists
    #hist, bins, patches = axh.hist(mag2, nbins, histtype='step', orientation='horizontal', log=True, color='black')
    #s_hist, s_bins, s_patches = axh.hist(s_mag2, s_nbins, histtype='step', orientation='horizontal', log=True, color='red')
    #axh.semilogx(hist, bins[:-1], drawstyle='steps', color='black')
    bins = b[:-1]
    axh.semilogx(mhist, bins, drawstyle='steps', color='darkgreen',
                 alpha=0.4)
    axh.semilogx(chist, bins, drawstyle='steps', color='darkblue',
                 alpha=0.4)

    kwargs['color'] = 'black'
    axh.annotate(r'$p=%.3f$' % p_value, xy=(.9, .95),
                 xycoords='axes fraction',
                 ha='right', **kwargs)
    '''
    # no formatters on mid and right plots
    [ax.yaxis.set_major_formatter(NullFormatter()) for ax in (axh, axm)]

    filename = '_'.join((gal.target, sgal.model_name, gal.filter1, gal.filter2))
    #plt.savefig(os.path.join(plt_dir, filename+'_LF.ps'))
    plt.savefig(os.path.join(plt_dir, sgal.mix, sgal.model_name, filename+'_LF.png'))
    print 'Wrote ', filename+'_LF.png'
    plt.close()
    return filename+'_LF.png'

if __name__ == '__main__':
    print 'use galaxy_test.py'
