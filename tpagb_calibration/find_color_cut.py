"""
The idea here was to stack all the data and find the color cuts well... that only
works with RGB of the same Z. doh.
"""
import os

import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import ResolvedStellarPops as rsp
from ResolvedStellarPops.tpagb_path_config import tpagb_path

data_loc = os.path.join(tpagb_path, 'SNAP/data/angst_no_trim/')
abs_mags = rsp.fileio.get_files(data_loc, '*dat')

wfpc2s = [a for a in abs_mags if 'wfpc2' in a]
acss =  [a for a in abs_mags if 'acs' in a]
instru = ['acs', 'wfpc2']

for i, cam in enumerate([acss, wfpc2s]):
    filter1s = np.array([os.path.split(a)[1].split('_')[5].split('-')[0]
                         for a in cam])
    filts = np.unique(filter1s)

    for filt in filts:
        gals = [a for a in cam if filt in a]
        mag1, mag2, Mag1, Mag2 = np.concatenate([np.genfromtxt(gal)
                                                 for gal in gals]).T
        fig, ax = plt.subplots()

        hist2d(Mag1 - Mag2, Mag2, bins=[100,100], cmap=plt.cm.Blues,
                       norm=mpl.colors.LogNorm())
        colorbar()
        ax.set_ylim(ax.get_ylim()[::-1])
        ax.set_ylabel(r'$F814W$', fontsize=20)
        ax.set_xlabel(r'${}-F814W$'.format(filt).upper(), fontsize=20)
        plt.tick_params(labelsize=16)
        ax.set_title('${}$'.format(instru[i].upper()))
        plt.savefig('alldata_{1}_{0}-f814w_hess.png'.format(filt, instru))
        
    
    
F606W - F814W wfpc2:
mag_bright = -2
mag_faint = -1
col_min = 0.5
col_max = 0.9

F435W - F814W acs:
mag_bright = -3
mag_faint = -2
col_min = 2.16
col_max = 2.3

F475W - F814W acs:
mag_bright = -3
mag_faint = -2
col_min = 2.16
col_max = 2.3