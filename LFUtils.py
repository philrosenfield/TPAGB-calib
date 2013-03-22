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


def plot_lines(axs, xrange, yval):
    [ax.plot((xrange), (yval, yval), color='black') for ax in axs]
    return


def plot_numbs(ax, item, xpos, ypos, **kwargs):
    ax.annotate(r'$%i$' % item, xy=(xpos, ypos), ha='left',
                **kwargs)
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


if __name__ == '__main__':
    print 'use galaxy_test.py'
