#!/usr/bin/env python
import argparse
import matplotlib as mpl

import matplotlib.pylab as plt
import numpy as np
import os
import sys
import time

import ResolvedStellarPops as rsp
from astropy.io import fits
from ResolvedStellarPops.tpagb_path_config import tpagb_path

from ..TPAGBparams import snap_src

data_loc = os.path.join(snap_src, 'data', 'galaxies')
    
def move_on(ok, msg='0 to move on: '):
    ok = int(raw_input(msg))
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
            ok = move_on(ok)

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
    here = os.getcwd()
    os.chdir(target)
    # read match file for mag1, mag2
    phot = glob.glob1('.', '*match')
    params = glob.glob1('.', '*param')
    for i, param in enumerate(params):
        # read param file
        mag1, mag2 = np.genfromtxt(phot[i], unpack=True)
        col = mag1 - mag2

        lines = open(param, 'r').readlines()
        colmin, colmax = map(float, lines[4].split()[3:-1])
        mag1min, mag1max = map(float, lines[5].split()[:-1])
        #mag2min, mag2max = map(float, lines[5].split()[:-1])
        # click around
        fig, ax = plt.subplots()
        ax.plot(col, mag2, ',', color='k', alpha=0.2)
        ax.set_ylim(mag1max, mag1min)
        ax.set_xlim(colmin, colmax)
    
        ok = 1
        while ok == 1:
            print 'click '
            pts = np.asarray(plt.ginput(n=4, timeout=-1))
            exclude_gate = '1 {} 0 \n'.format(' '.join(['%.4f' % p for p in pts.flatten()]))
            pts = np.append(pts, pts[0]).reshape(5,2)
            ax.plot(pts[:,0], pts[:,1], color='r', lw=3, alpha=0.3)
            plt.draw()
            ok = move_on(0)
        lines[7] = exclude_gate
        # not so simple ... need them to be parallelograms.
        # PASS!
        
        # write new param file with exclude/include gate
        os.system('mv {0} {0}_bkup'.format(param))
        with open(param, 'w') as outp:
            [outp.write(l) for l in lines]
        print('wrote %s' % param)
    
    os.chdir(here)


def laad_obs(target, optfilter1=''):
    """load in NIR and OPT galaxy"""
    nirgalname, = rsp.fileio.get_files(data_loc, '*{}*fits'.format(target.upper()))
    optgalname, = rsp.fileio.get_files(data_loc, ('*{}*{}*fits'.format(target, optfilter1).lower()))

    nirphot = fits.getdata(nirgalname)

    optphot = fits.getdata(optgalname)
    return optphot, nirphot

def find_normalization_limits(target, optfilter1=''):
    optphot, nirphot = laad_obs(target, optfilter1=optfilter1)
    phot_ext = 'acs'
    if 'wfpc2' in optphot:
        phot_ext = 'wfpc2'
    
    optcolmin, optcolmax = find_match_limits(optphot, phot_ext, color_only=True)
    nircolmin, nircolmax = find_match_limits(nirphot, 'ir', color_only=True)
    print ','.join('%.2f' % m for m in [optcolmin, optcolmax, nircolmin, nircolmax])

def match_limits(color_only=False, data_file='snap_galaxies.dat',
                 target=None):
    plt.ion()
    new_lines = '# target comp_nir1 comp_nir2 comp_opt1 comp_opt2 Av mTRGB mTRGBerr dmod colmin colmax mag1max mag2max opt1 opt2 opt_phot opt_fake\n'
    data_loc = os.path.join(tpagb_path, 'SNAP/data/angst_no_trim')
    lines = open(os.path.join(tpagb_path, 'SNAP/tables/{}'.format(data_file)), 'r').readlines()

    if target is not None:
        lines = [l for l in lines if target in l]

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

    parser.add_argument('-t', '--target', type=str, default=None,
                        help='name of one target in the data_file')

    parser.add_argument('-e', '--exgates', action='store_true',
                        help='Find exclude gates instead')

    parser.add_argument('-n', '--norm', action='store_true',
                        help='Find normalization color limits')
    
    parser.add_argument('-f', '--optfilter1', type=str, default='',
                        help='optical V filter')
    
    args = parser.parse_args(sys.argv[1:])
    color_only = args.color_only
    if args.pdb:
        import pdb
        pdb.set_trace()
    if args.exgates:
        assert args.target is not None, \
            'Must supply target if finding exclude gates'
        find_gates(args.target)
    elif args.norm:
        find_normalization_limits(args.target, optfilter1='')
    else:
        match_limits(color_only=color_only, data_file=args.data_file, target=args.target)
