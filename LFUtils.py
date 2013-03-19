#
#  LumnFncfromCMD.py
#
#  Created by Philip Rosenfield on 11/22/10.
#
import os
import numpy as np
import brewer2mpl
import matplotlib.pyplot as plt
from matplotlib.ticker import NullFormatter
from matplotlib.patheffects import withStroke
from matplotlib import rc
myeffect = withStroke(foreground="w", linewidth=3)
kwargs = dict(path_effects=[myeffect])
rc('text', usetex=True)
#from ResolvedStellarPops.math_utils import bayesian_blocks



from scipy.stats import ks_2samp



# Phil's functions:
import GenUtils
from TrilegalUtils import get_label_stage
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




def brighter(mag2, trgb, inds=None):
    ''' number of stars brighter than trgb, make sure mag2 is
        the same filter as trgb!'''
    i, = np.nonzero(mag2 < trgb)
    if inds is not None:
        i = np.intersect1d(i, inds)
    return i




def get_chi2(obs, exp):
    return np.sum((obs - exp) ** 2 / exp)




def between(arr, mdim, mbrt, inds=None):
    i, = np.nonzero((arr < mdim) & (arr > mbrt))
    if inds is not None:
        i = np.intersect1d(i, inds)
    return i




def calc_LF(gal, sgal, maglims, res=0.1, normalize=True):
    '''
    usage:
    gal: galaxy object
    sgal: simgalaxy object
    maglims: tuple (region in mag to normalize)
    res: LF bin resolution [0.1]
    normalize: bool [True] do the normalization



    adds attributes to gal and sgal:
        irgb: indices of "rgb" stars (and agb stars, relys on tagged data)
        iagb: indices of "rgb" stars brighter than trgb
        bins: LF bins
        LF: Luminosity function



    ------
    adds attributes to gal:
    rgb_norm: indices of rgb stars between maglims



    ------
    adds attributes to sgal:



    norm:
        (only if normalization)
        indices of all stars in same cmd space as data's rgb_norm



    rel_ind:
        if normalize = True: indices of simulation cmd that are picked so they
                             mimic the number of rgb stars between maglims
        else: the entire array.



    normalization:
        (only if normalization)
        scaling used to determine ind and scale LFs (added attributes to
        objects)



    LFn:
        (only if normalization) LF scaled by normalization
    rel_rgb:
        indices of rgb stars randomly picked from normalization
    rel_agb:
        indices of agb stars randomly picked from normalization



    returns:
    p_value: ks test probability comparing simulated stars brighter than trgb
             to data of the same



    '''



    # data
    mag = gal.mag2
    color = gal.color
    dim, bright = maglims
    # when hand picking cmd_regions, RGB includes AGB!
    gal.irgb = gal.stage_inds('RGB')
    gal.iagb = brighter(mag, gal.trgb, inds=gal.irgb)



    # simulation
    smag = sgal.ast_mag2
    scolor = sgal.ast_color
    sgal.irgb = sgal.stage_inds('RGB')



    # LF resolution in dex
    nbins = (mag.max() - mag.min()) / res



    # the LFs
    gal.LF, gal.bins = np.histogram(mag, nbins)
    sgal.LF, sgal.bins = np.histogram(smag, bins=gal.bins)



    if normalize is True:
        # The RGB stars used for normalization
        gal.rgb_norm = between(mag, dim, bright, inds=gal.irgb)



        # The Simulated RGB stars
        #!! Not used!!
        #sgal.rgb_in_norm = between(smag, dim, bright, inds=sgal.irgb)



        # all sim stars in cmd region of rgb data norm stars.
        sgal.norm = GenUtils.inside(color[gal.rgb_norm], mag[gal.rgb_norm],
                                    scolor, smag)



        # here is the normalization!!
        sgal.normalization = float(gal.rgb_norm.size)/float(sgal.norm.size)
        if sgal.normalization > 0.75:
            print sgal.normalization
            print 'Too few model stars'



        # normalize the simulated LF
        sgal.LFn = sgal.LF*sgal.normalization



        # for plotting random simulated stars to match number of obs.
        rands = np.random.random(smag.size)
        ind, = np.nonzero(rands < sgal.normalization)
        sgal.rel_ind = ind
    else:
        ind = range(smag.size)



    # simulated normalized RGB stars
    sgal.rel_rgb = between(smag[ind], dim, gal.trgb)



    # simulated normalized stars brighter than trgb.
    sgal.rel_agb = brighter(smag[ind], gal.trgb - sgal.count_offset)



    KS_D, p_value = ks_2samp(mag[gal.iagb], smag[sgal.rel_agb])



    # want to add some lines for the flux and the mass loss rates.
    return p_value




def setup_lfplot(gal):
    # plot limits
    left, width = 0.1, 0.312
    bottom, height = 0.1, 0.8
    left_m = left + width + 0.01  # for model
    left_h = left + 2 * width + 0.02  # for hist (width is set hist_axis)



    # plot and fig sizes
    fig = plt.figure(1, figsize=(8, 8))



    data_axis = [left, bottom, width, height]
    model_axis = [left_m, bottom, width, height]
    hist_axis = [left_h, bottom, 0.2, height]



    axData = plt.axes(data_axis)
    axModel = plt.axes(model_axis)
    axHist = plt.axes(hist_axis)



    # titles
    axData.set_title(r'$\rm{data}$', color='black')
    axModel.set_title(r'$\rm{model}$', color='red')
    axData.set_xlabel(r'$\rm{%s-%s}$' % (gal.filter1, gal.filter2), size=20)
    axData.set_ylabel(r'$\rm{%s}$' % gal.filter2, size=20)
    axModel.set_xlabel(axData.get_xlabel(), size=20)
    axHist.set_xlabel(r'$\#$', size=20)



    for ax in [axData, axModel, axHist]:
        ax.set_xlim((-0.5, gal.color.max()))  # model and data x limits here
        ax.set_ylim((25.5, 18.5))  # set all y limits here
    axHist.set_xlim((1, 50000))



    return fig, axData, axModel, axHist




def plot_lines(axs, xrange, yval):
    [ax.plot((xrange), (yval, yval), color='black') for ax in axs]
    return




def plot_numbs(ax, item, xpos, ypos, **kwargs):
    ax.annotate(r'$%i$' % item, xy=(xpos, ypos), ha='left', size=20,
                **kwargs)
    return




def diagnostic_cmd(sgal, trgb, figname=None, inds=None):
    if inds is not None:
        ustage = np.unique(sgal.stage[inds])
    else:
        ustage = np.unique(sgal.stage)
    nplots = ustage.size+1.
    cols = brewer2mpl.get_map('Paired', 'qualitative', len(ustage)).mpl_colors
    subplots_kwargs = {'sharex': 1, 'sharey': 1, 'figsize': (12, 8)}
    j = 0
    for color, mag2 in zip((sgal.color, sgal.ast_color),
                           (sgal.mag2, sgal.ast_mag2)):
        if inds is not None:
            stage = sgal.stage[inds]
        fig, (axs) = setup_multiplot(nplots, **subplots_kwargs)



        for ax in axs.ravel():
            ax.set_xlim(-0.5, sgal.color.max())
            ax.set_ylim(25.5, 18.5)



        ax0, cols = colorplot_by_stage(axs.ravel()[0], color, mag2, '.', stage,
                                       cols=cols)
        i = 0
        for ax, st in zip(axs.ravel()[1:], ustage):
            plot_lines([ax], ax.get_xlim(), trgb)
            label = get_label_stage(int(st))
            ind = sgal.stage_inds(label)
            if inds is not None:
                ind = list(set(ind) & set(inds))
            if len(ind) == 0:
                continue
            ax.plot(color[ind], mag2[ind], '.', color=cols[i], mew=0,
                    label='N=%i' % len(ind))
            kwargs['color'] = 'black'
            text_offset = 0.02
            xpos = ax.get_xlim()[0]+2*text_offset
            plot_numbs(ax, brighter(mag2, trgb-sgal.count_offset, inds=ind).size, xpos,
                       trgb-text_offset, **kwargs)
            ax.set_title(label, **{'color': cols[i]})
            i += 1
            ax.legend(loc=1, numpoints=1, frameon=False)
        if figname:
            if j == 0:
                extra = ''
            else:
                extra = '_spread'
            plt.savefig(figname.replace('.png', '%s.png' % extra))
            print 'wrote %s' % figname.replace('.png', '%s.png' % extra)
            plt.close()
        else:
            plt.show()
        j += 1
    return figname.replace('.png', '%s.png' % extra)




def setup_multiplot(nplots, **subplots_kwargs):
    nx = np.round(np.sqrt(nplots))
    nextra = nplots-nx**2
    ny = nx
    if nextra > 0:
        ny += 1
    nx = int(nx)
    ny = int(ny)



    fig, axs = plt.subplots(nx, ny, **subplots_kwargs)



    return fig, axs




def colorplot_by_stage(ax, x, y, marker, stages, cols=None):
    # inds from calc_LFIR are based on only resolved stars.



    if cols is None:
        cols = discrete_colors(len(np.unique(stages)))
    for i, s in enumerate(np.unique(stages)):
        ind, = np.nonzero(stages == s)
        if ind.size == 0:
            continue
        ax.plot(x[ind], y[ind], marker, color=cols[i], mew=0)
    return ax, cols




def make_title(gal, fig):
    text_kwargs = {'ha': 'center', 'va': 'top', 'size': 20}
    if np.isfinite(gal.z):
        title = r'$\rm{%s\ m}_{TRGB}=%.3f\ Z=%.4f$' % (gal.target, gal.trgb,
                                                       gal.z)
    else:
        title = r'$\rm{%s\ m}_{TRGB}=%.3f\ Z=...$' % (gal.target, gal.trgb)



    fig.text(0.5, 0.96, title, **text_kwargs)
    return




def plot_LFIR(gal, sgal, p_value, maglims):
    dim, bright = maglims



    mstars = list(set(sgal.imstar) & set(sgal.rel_ind))
    cstars = list(set(sgal.icstar) & set(sgal.rel_ind))



    #nbright_rgb = bright_rgb.size
    bright_mstars = brighter(sgal.ast_mag2, gal.trgb-sgal.count_offset,
                             inds=mstars)
    bright_cstars = brighter(sgal.ast_mag2, gal.trgb-sgal.count_offset,
                             inds=cstars)



    nbright_rgb = sgal.rel_agb.size - len(bright_mstars) - len(bright_cstars)



    # model used for normalization
    norm_inds = list(set(sgal.norm) & set(sgal.rel_ind))



    if sgal.rel_agb.size - len(bright_mstars) - len(bright_cstars) - nbright_rgb != 0.:
        print ''
        print sgal.rel_agb.size, len(bright_mstars) - len(bright_cstars) - nbright_rgb
        print ''



    #nRGBs = len(sgal.ast_color[sgal.rel_agb]) - len(mstars) - len(cstars)



    mhist, b = np.histogram(sgal.ast_mag2[mstars], bins=sgal.bins)
    chist, b = np.histogram(sgal.ast_mag2[cstars], bins=sgal.bins)



    fig, axData, axModel, axHist = setup_lfplot(gal)



    make_title(gal, fig)



    # plot data
    axData.plot(gal.color, gal.mag2, '.', mew=0, color='grey', mec='grey')
    # data used for normalization
    axData.plot(gal.color[gal.rgb_norm], gal.mag2[gal.rgb_norm], '.', mew=0,
                color='black')



    # plot model
    axModel.plot(sgal.ast_color[sgal.rel_ind], sgal.ast_mag2[sgal.rel_ind],
                 '.', color='red', mew=0, label=r'Other$=%i$' % (nbright_rgb))



    axModel.plot(sgal.ast_color[norm_inds], sgal.ast_mag2[norm_inds], '.',
                 color='black', mew=0)  # , label=r'$RGB=%i$'%(len(norm_inds)))



    axModel.plot(sgal.ast_color[mstars], sgal.ast_mag2[mstars], '.',
                 color='darkblue', mew=0, label=r'$M=%i$' % (len(bright_mstars)))
    #             color='darkgreen', mew=0, label=r'$M=%i$'%(len(mstars)))



    axModel.plot(sgal.ast_color[cstars], sgal.ast_mag2[cstars], '.',
                 color='darkblue', mew=0, label=r'$C=%i$' % (len(bright_cstars)))
    #             color='darkblue', mew=0, label=r'$C=%i$'%(len(cstars)))
    axModel.legend(frameon=False, loc=2, numpoints=1)



    #[plot_lines([axData, axModel], axData.get_xlim(), m) for m in (dim, bright)]
    # also?
    plot_lines([axData, axModel], axData.get_xlim(), gal.trgb)
    plot_lines([axData, axModel], axData.get_xlim(), gal.trgb-sgal.count_offset)



    kwargs['color'] = 'black'
    text_offset = 0.02
    xpos = axData.get_xlim()[0] + 2 * text_offset
    yposs = np.asarray(maglims) - text_offset
    ypos = gal.trgb-text_offset
    plot_numbs(axData, gal.iagb.size-sgal.count_offset, xpos, ypos, **kwargs)
    plot_numbs(axData, gal.rgb_norm.size, xpos, yposs[0], **kwargs)



    kwargs['color'] = 'red'
    xpos = axModel.get_xlim()[0] + text_offset
    plot_numbs(axModel, sgal.rel_agb.size, xpos, ypos, **kwargs)
    plot_numbs(axModel, len(norm_inds), xpos, yposs[0], **kwargs)



    # plot hists
    #hist, bins, patches = axHist.hist(mag2, nbins, histtype='step', orientation='horizontal', log=True, color='black')
    #s_hist, s_bins, s_patches = axHist.hist(s_mag2, s_nbins, histtype='step', orientation='horizontal', log=True, color='red')
    #axHist.semilogx(hist, bins[:-1], drawstyle='steps', color='black')
    bins = b[:-1]
    axHist.semilogx(mhist, bins, drawstyle='steps', color='darkgreen',
                    alpha=0.4)
    axHist.semilogx(chist, bins, drawstyle='steps', color='darkblue',
                    alpha=0.4)
    axHist.semilogx(gal.LF, bins, drawstyle='steps', color='black', lw=2)
    axHist.semilogx(sgal.LFn, bins, drawstyle='steps', color='red', lw=2)
    kwargs['color'] = 'black'
    axHist.annotate(r'$p=%.3f$' % p_value, xy=(.9, .95),
                    xycoords='axes fraction',
                    ha='right', size=20, **kwargs)



    # no formatters on mid and right plots
    nullfmt = NullFormatter()  # no labels
    [ax.yaxis.set_major_formatter(nullfmt) for ax in (axHist, axModel)]



    filename = '_'.join((gal.target, sgal.model_name, gal.filter1, gal.filter2))
    #plt.savefig(os.path.join(plt_dir, filename+'_LF.ps'))
    plt.savefig(os.path.join(plt_dir, sgal.mix, sgal.model_name, filename+'_LF.png'))
    print 'Wrote ', filename+'_LF.png'
    plt.close()
    return filename+'_LF.png'



if __name__ == '__main__':
    print 'use galaxy_test.py'
