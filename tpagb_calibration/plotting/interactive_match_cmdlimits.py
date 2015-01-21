#!/usr/bin/python
import argparse
import matplotlib.pylab as plt
import numpy as np
import os
import sys
import time

from astropy.io import fits
from ResolvedStellarPops.tpagb_path_config import tpagb_path

def move_on(ok, msg='0 to move on: '):
    #ok = int(raw_input(msg))
    time.sleep(1)
    return ok

def find_match_limits(phot, phot_ext, comp1=99., comp2=None, color_only=False,
                      xlim=None, ylim=None):
    """
    click color limits on a cmd and mag1 mag2 limits on a plot of mag1 vs mag2
    """
    mag1 = phot['mag1_%s' % phot_ext]
    mag2 = phot['mag2_%s' % phot_ext]
    col = mag1 - mag2

    fig, ax = plt.subplots()
    ax.plot(col, mag2, 'o', color='k', ms=3, alpha=0.3, mec='none')
    if xlim is not None:
        ax.set_xlim(xlim)
    if ylim is not None:
        ax.set_ylim(ylim)
    else:
        ax.set_ylim(ax.get_ylim()[::-1])
    
    if comp2 is not None:
        ax.hlines(comp2, *ax.get_xlim())
    else:
        comp2 = 99.

    ok = 1
    while ok == 1:
        print 'click on color min then color max'
        pts = plt.ginput(2, timeout=-1)
        colmin, colmax = [pts[i][0] for i in range(2)]
        ax.vlines(colmin, *ax.get_ylim())
        ax.vlines(colmax, *ax.get_ylim())
        plt.draw()
        ok = move_on(0)

    plt.close()

    inds, = np.nonzero((col < colmax) & (col > colmin))
    data = (colmin, colmax)
    if not color_only:
        fig, ax = plt.subplots()
        ax.plot(mag1, mag2, '.', color='k')
        ok = 1
        while ok == 1:
            print 'click the bright mag value of mag1 and mag2, click a second time to finish'
            pts = plt.ginput(2, timeout=-1)
            mag1max, mag2max = pts[0]
            ax.plot(mag1max, mag2max, 'o', color='r')
            plt.draw()
            ok = move_on(0)

        plt.close()

        inds, = np.nonzero((mag1 < comp1) & (mag1 > mag1max) &
                           (mag2 < comp2) & (mag2 > mag2max) &
                           (col < colmax) & (col > colmin))

    fig, ax = plt.subplots()
    ax.plot(col, mag2, '.', color='k')
    ax.plot(col[inds], mag2[inds], '.', color='r')
    ax.set_ylim(ax.get_ylim()[::-1])
    ax.hlines(comp2, *ax.get_xlim(), lw=2)
    ax.vlines(colmin, *ax.get_ylim(), lw=2)
    ax.vlines(colmax, *ax.get_ylim(), lw=2)
    if not color_only:
        ax.hlines(mag2max, *ax.get_xlim(), lw=2)
        data = (colmin, colmax, mag1max, mag2max)
    
    plt.draw()

    return data

def find_gates(target):
    import glob
    target = 'ddo82'
    here = os.getcwd()
    os.chdir(target)
    # read match file for mag1, mag2
    phot, = glob.glob1('.', '*match')
    mag1, mag2 = np.genfromtxt(phot, unpack=True)
    col = mag1 - mag2
    # read param file
    param, = glob.glob1('.', '*param')
    lines = open(param, 'r').readlines()
    colmin, colmax = map(float, lines[4].split()[3:-1])
    mag1min, mag1max = map(float, lines[5].split()[:-1])
    #mag2min, mag2max = map(float, lines[5].split()[:-1])
    # click around
    fig, ax = plt.subplots()
    ax.plot(col, mag2, ',', color='k', alpha=0.1)
    ax.set_ylim(mag1max, mag1min)
    ax.set_xlim(colmin, colmax)

    ok = 1
    while ok == 1:
        print 'click '
        pts = plt.ginput(-1, timeout=-1)
        mag1max, mag2max = pts[0]
        ax.plot(mag1max, mag2max, 'o', color='r')
        plt.draw()
        ok = move_on(0)
    # not so simple ... need them to be parallelograms.

    # write new param file with exclude/include gate
    os.chdir(here)
    

def match_limits(color_only=False, data_file='snap_galaxies.dat'):
    plt.ion()
    new_lines = '# target comp_nir1 comp_nir2 comp_opt1 comp_opt2 Av mTRGB mTRGBerr dmod colmin colmax mag1max mag2max opt1 opt2 opt_phot opt_fake\n'
    data_loc = os.path.join(tpagb_path, 'SNAP/data/angst_no_trim')
    lines = open(os.path.join(tpagb_path, 'SNAP/tables/{}'.format(data_file)), 'r').readlines()
    for line in lines:
        if line.startswith('#'):
            continue
        target, comp_nir1, comp_nir2, comp_opt1, comp_opt2, Av, mTRGB, mTRGBerr, dmod, cmin, cmax, opt1, opt2, filter1, filter2, opt_phot, opt_fake = line.split()
        opt_phot = os.path.join(data_loc, opt_phot)
        opt_fake = os.path.join(data_loc, opt_fake)
        assert os.path.isfile(opt_phot), 'no opt phot'
        assert os.path.isfile(opt_fake), 'no opt fake'
        print target
        phot = fits.getdata(opt_phot)
        phot_ext = 'acs'
        if 'wfpc2' in opt_phot:
            phot_ext = 'wfpc2'

        ok = 1
        while ok == 1:
            data = find_match_limits(phot, phot_ext, comp1=float(comp_opt1),
                                     comp2=float(comp_opt2),
                                     xlim=(float(cmin), float(cmax)),
                                     color_only=color_only)
            ok = move_on(0)

        plt.close()

        partial_line = ' '.join([target, comp_nir1, comp_nir2, comp_opt1,
                                 comp_opt2, Av, mTRGB, mTRGBerr, dmod])
        end_line = ' '.join([opt1, opt2, opt_phot, opt_fake])

        if color_only:
            colmin, colmax = data
            data_str = '%.2f %.2f' % (colmin, colmax)
        else:
            colmin, colmax, mag1max, mag2max = data
            data_str = '%.2f %.2f %.2f %.2f' % (colmin, colmax, mag1max, mag2max)
        
        new_line = '%s %s %s \n' % (partial_line, data_str, end_line)

        print new_line
        new_lines += new_line

    print new_lines
    
if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Find color mag limits of CMDs interactively")

    parser.add_argument('-c', '--color_only', action='store_true',
                        help='skip the magnitude finding')

    parser.add_argument('-d', '--data_file', type=str, default='snap_galaxies.dat',
                        help='data table in [tpagb_path]/SNAP/tables to use')

    parser.add_argument('-v', '--pdb', action='store_true',
                        help='toggle debugging')

    args = parser.parse_args(sys.argv[1:])
    color_only = args.color_only
    if args.pdb:
        import pdb
        pdb.set_trace()
    match_limits(color_only=color_only, data_file=args.data_file)
