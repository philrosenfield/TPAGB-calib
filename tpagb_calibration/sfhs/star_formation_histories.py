from __future__ import print_function
import logging
import os

import matplotlib.pylab as plt
import numpy as np
import ResolvedStellarPops as rsp

from ResolvedStellarPops.convertz import convertz

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

__all__ = ['StarFormationHistories', 'parse_sfh_data']

def parse_sfh_data(filename, file_origin, frac=0.2):
    '''
    parse match sfh into a np.recarray
    ARGS:
    filename: a match sfh file to parse, needs to have at least:
              dtype = [('lagei', '<f8'),
                       ('lagef', '<f8'),
                       ('sfr', '<f8'),
                       ('sfr_errp', '<f8'),
                       ('sfr_errm', '<f8'),
                       ('mh', '<f8'),
                       ('mh_errp', '<f8'),
                       ('mh_errm', '<f8'),
                       ('mh_disp', '<f8')]
    file_origin: 'match' or 'match-grid'
        if 'match': reads the binned sfh as takes sfr_errs as is
        if 'match-grid': will overwrite sfr_errs.
            if sfr = 0: sfr_errp = min(sfr) * frac
            if sfr != 0: sfr_errp = sfr_errm = sfr * frac
    frac: multiplicitve factor to set sfr_err, only used if file_origin
        is set to match-grid.
    RETURNS:
    np.recarray of the sfh file

    use file_origin='match' only if Hybrid MonteCarlo has been run.
    '''
    if 'match-grid' == file_origin.lower() or 'match-hmc' == file_origin.lower():
        data = rsp.match.utils.read_binned_sfh(filename)
    elif 'match-old' == file_origin.lower():
        data = rsp.match.utils.read_match_old(filename)
    else:
        logger.error('please add a new data reader')

    if 'grid' in file_origin.lower():
        logger.warning('CONSTANT %f SFR ERROR' % frac)
        # at least have a uniform error that is 20% the smallest sfr.
        data.sfr_errp = np.min(data.sfr[data.sfr > 0]) * frac
        # now take 10% as nominal error in each sfr bin.
        data.sfr_errp[data.sfr > 0] = data.sfr[data.sfr > 0] * frac
        data.sfr_errm[data.sfr > 0] = data.sfr[data.sfr > 0] * frac
    return data


class StarFormationHistories(object):
    '''Make TRILEGAL star formation history files from MATCH'''
    def __init__(self, sfh_file, file_origin, sfr_files=None, sfr_file_loc=None,
                 sfr_file_search_fmt=None, sfh_ext='.sfh'):
        self.base, self.name = os.path.split(sfh_file)
        self.data = parse_sfh_data(sfh_file, file_origin)
        self.file_origin = file_origin
        self.sfr_files = sfr_files
        self.sfh_ext = sfh_ext
        if sfr_file_loc is not None:
            if sfr_file_search_fmt is not None:
                self.sfr_files = rsp.fileIO.get_files(sfr_file_loc,
                                                      sfr_file_search_fmt)

    def random_draw_within_uncertainty(self, attr, npoints=2e5):
        '''
        randomly draw values within uncertainty for an array

        ARGS:
        attr: string name of the array that also has attr_errm and
              attr_errp (p and m are important due to the sign).
              attr_errm: - err associated with each point on attr array
              attr_errp: same as attr_errm but the + err

        npoints: number of points to populate gaussian to sample

        RETURNS:
        array of values randomly picked within the uncertainties

        If errm and errp are equal, just returns a randomly chosen
        point (of npoints) of a gaussian with mean attr and
        sigma=attr_errm

        If not, will stick to gaussians together at attr using
        sigma=attr_errm and sigm=attr_errp and returning a random value
        from there.

        If one of the err values is zero, will just use the other half
        of the gaussian.

        If they are both zero, well, just returns attr.
        '''
        if attr == 'mh':
            logger.warning('this method was designed for sfr, not [M/H]')
        # load in values this way in case I want to move this to its own
        # function
        if hasattr(self.data, attr):
            val_arr = self.data.__getattribute__(attr)
            errm_arr = self.data.__getattribute__('%s_errm' % attr)
            errp_arr = self.data.__getattribute__('%s_errp' % attr)
        else:
            val_arr = self.__getattribute__(attr)
            errm_arr = self.__getattribute__('%s_errm' % attr)
            errp_arr = self.__getattribute__('%s_errp' % attr)
        rand_arr = np.array([])
        # don't want negative sfr values. If not sfr, don't care.
        if attr == 'sfr':
            lowlim = 0
        else:
            lowlim = -np.inf

        for val, errm, errp in zip(val_arr, errm_arr, errp_arr):
            if errp == errm and errp > 0:
                # even uncertainties, easy.
                new_arr = np.random.normal(val, errp, npoints)
            elif errp != 0 and errm != 0:
                # stitch two gaussians together
                pos_gauss = np.random.normal(val, errp, npoints)
                neg_gauss = np.random.normal(val, errm, npoints)
                new_arr = np.concatenate([pos_gauss[pos_gauss >= val],
                                          neg_gauss[neg_gauss <= val]])

            elif errp == 0 and errm != 0:
                # no positive uncertainties
                neg_gauss = np.random.normal(val, errm, npoints)
                new_arr = neg_gauss[neg_gauss <= val]
            elif errp != 0 and errm == 0:
                # no negative uncertainties
                pos_gauss = np.random.normal(val, errp, npoints)
                new_arr = pos_gauss[pos_gauss >= val]
            else:
                # um.. no errors, why was this called
                logger.warning('no uncertainties')
                new_arr = np.ones(4) * val
            new_arr = new_arr[new_arr > lowlim]
            rand_arr = np.append(rand_arr, np.random.choice(new_arr))
        return rand_arr

    def interp_null_values(self):
        '''
        If there is no SF, there is still some +err in SF. However, M/H is
        not constrained so is set to 0. Here we fill in values of M/H by
        interpolating the entire M/H vs age, this should be used as mean value
        with vdisp to the be the sigma in the gaussian distribution.

        I think it's reasonable since this is really just finding the -zinc
        law that MATCH assumes.
        '''
        somesf, = np.nonzero(self.data.sfr != 0)
        #ff = interp1d(self.data.lagei[somesf], self.data.mh[somesf],
        #              bounds_error=False)
        _, mh_interp = rsp.math_utils.extrap1d(self.data.lagei[somesf],
                                               self.data.mh[somesf],
                                               self.data.lagei)

        self.mh_interp = mh_interp
        return mh_interp

    def make_trilegal_sfh(self, random_sfr=False, random_z=False,
                          zdisp=True, outfile='default', dry_run=False):
        '''
        turn binned sfh in to trilegal sfh
        random_sfr:
            calls random_draw_within_uncertainty
        random_z:
        '''
        # In MATCH [M/H] = log(Z/Zsun) with Zsun = 0.02 (see MATCH's makemod.cpp)
        # It doesn't matter if this is "correct". Stellar models have absolute Z.
        # Zsun is just a scaling that needs to be undone from MATCH to here.
        zsun = 0.02

        if self.file_origin == 'match-old':
            random_sfr = False
            random_z = False
            zdisp = False
        if outfile == 'default':
            outfile = os.path.join(self.base,
                                   self.name.replace(self.sfh_ext, '.tri.dat'))
        if dry_run is True:
            return outfile

        age1a = 10 ** (self.data.lagei)
        age1p = 1.0 * 10 ** (self.data.lagei + 0.0001)
        age2a = 1.0 * 10 ** self.data.lagef
        age2p = 1.0 * 10 ** (self.data.lagef + 0.0001)

        if random_sfr is False:
            sfr = self.data.sfr
        else:
            sfr = self.random_draw_within_uncertainty('sfr')


        if random_z is False:
            mh = self.data.mh
        else:
            # HACK. Not using mh errs from MATCH. Untrustworthy.
            # Shifting instead from within dispersion.
            self.interp_null_values()
            disp = np.median(self.data.mh_disp[np.nonzero(self.data.mh_disp)])/2.
            mh = self.mh_interp + np.random.normal(0, disp)

        metalicity = zsun * 10 ** mh

        if zdisp is True:
            zdisp = metalicity * np.median(self.data.mh_disp[np.nonzero(self.data.mh_disp)])
            #zdisp = self.data.mh_disp
            fmt = '%.4e %.3e %.4f %.4f \n'
        else:
            zdisp = [''] * len(mh)
            fmt = '%.4e %.3e %.4f %s\n'

        with open(outfile, 'w') as out:
            for i in range(len(sfr)):
                if sfr[i] == 0:
                    # this is just a waste of lines in TRILEGAL
                    continue
                out.write(fmt % (age1a[i], 0.0, metalicity[i], zdisp[i]))
                out.write(fmt % (age1p[i], sfr[i], metalicity[i], zdisp[i]))
                out.write(fmt % (age2a[i], sfr[i], metalicity[i], zdisp[i]))
                out.write(fmt % (age2p[i], 0.0, metalicity[i], zdisp[i]))
        return outfile

    def load_random_arrays(self, attr_str):
        if 'sfr' in attr_str:
            col = 1
        if 'mh' in attr_str or 'feh' in attr_str:
            col = 2
        val_arrs = [np.genfromtxt(s, usecols=(col))[1::4]
                    for s in self.sfr_files]
        if attr_str == 'mh':
            val_arrs = np.array([10**(val_arr/.2) for val_arr in val_arrs])
        if attr_str == 'feh':
            val_arrs = np.array([convertz.convertz(feh=feh)[1]
                                 for feh in val_arrs])

        return val_arrs

    def plot_sfh(self, attr_str, ax=None, outfile=None, yscale='linear',
                 plot_random_arrays_kw=None, errorbar_kw=None,
                 twoplots=False):
        '''
        plot the data from the sfh file.
        '''
        plot_random_arrays_kw = plot_random_arrays_kw or {}

        # set up errorbar plot
        errorbar_kw = errorbar_kw or {}
        errorbar_default = {'linestyle': 'steps-mid', 'lw': 3, 'color': 'darkred'}
        errorbar_kw = dict(errorbar_default.items() + errorbar_kw.items())

        # load the plotting values and their errors, this could be generalized
        # and passed ...
        val_arr = self.data.__getattribute__(attr_str)
        errm_arr = self.data.__getattribute__('%s_errm' % attr_str)
        errp_arr = self.data.__getattribute__('%s_errp' % attr_str)

        if 'sfr' in attr_str:
            ylab = '${\\rm SFR\ (10^3\ M_\odot/yr)}$'
            val_arr *= 1e3
            errm_arr *= 1e3
            errp_arr *= 1e3
            plot_random_arrays_kw['moffset'] = 1e3
        elif 'm' in attr_str:
            ylab = '${\\rm [M/H]}$'
        elif 'fe' in attr_str:
            ylab = '${\\rm [Fe/H]}$'

        # mask 0 values so there is a vertical line on the plot
        val_arr[val_arr==0] = 1e-15
        errm_arr[errm_arr==0] = 1e-15
        errp_arr[errp_arr==0] = 1e-15

        if twoplots is True:
            fig, axs = plt.subplots(figsize=(8, 8), nrows=2, sharex=True)
        else:
            fig, ax = plt.subplots(figsize=(8, 5))
            axs = [ax]
        for ax in axs:
            if not 'm' in attr_str:
                ax.errorbar(self.data.lagei, val_arr, [errm_arr, errp_arr],
                            zorder=100, **errorbar_kw)

            if len(plot_random_arrays_kw) > 0:
                # if loading the random arrays from files, need to give the
                # attribute to load.
                if plot_random_arrays_kw['from_files'] is True:
                    plot_random_arrays_kw['attr_str'] = attr_str
                self.plot_random_arrays(ax=ax, **plot_random_arrays_kw)
            ax.set_ylabel(ylab, fontsize=20)
            ax.set_xlim(7.9, 10.13)

        target = self.name.split('.')[0].upper().replace('-', '\!-\!')
        if twoplots is True:
            ax = axs[1]
            # lower plot limit doesn't need to be 1e-15...
            axs[0].set_ylim(1e-7, axs[0].get_ylim()[1])
            axs[0].set_yscale('log')
            fig.subplots_adjust(hspace=0.09)
        else:
            fig.subplots_adjust(bottom=0.15)

        ax.annotate('$%s$' % target, (0.02, 0.97), va='top',
                        xycoords='axes fraction', fontsize=16)
        ax.set_xlabel('$\log {\\rm Age (yr)}$', fontsize=20)

        plt.tick_params(labelsize=16)

        if outfile is not None:
            plt.savefig(outfile, dpi=150)
        return ax

    def plot_random_arrays(self, ax=None, val_arrs=None, from_files=False,
                           attr_str=None, yscale='linear', moffset=1.):
        '''
        val_arrs are random
        after making a bunch of arrays that sample the sfr or mh uncertainties
        plot up where they are.
        '''
        if from_files is True:
            val_arrs = self.load_random_arrays(attr_str)

        assert val_arrs is not None, 'either specify val_arrs or set from_files'

        if ax is None:
            fig, ax = plt.subplots(figsize=(12,12))

        [ax.errorbar(self.data.lagei, val_arrs[i] * moffset,
                     linestyle='steps-mid',
                     color='k', alpha=0.2) for i in range(len(val_arrs))]
        ax.set_yscale('linear')
        return ax

    def make_many_trilegal_sfhs(self, nsfhs=100, outfile_fmt='default',
                                random_sfr=True, random_z=False,
                                zdisp=True, dry_run=False):
        '''
        make nsfhs number of trilegal sfh input files.
        '''
        if outfile_fmt == 'default':
            outfile_fmt = self.name.replace(self.sfh_ext, '%03d.tri.dat')

        mk_tri_sfh_kw = {'random_sfr': random_sfr, 'random_z': random_z,
                         'zdisp': zdisp, 'dry_run': dry_run}

        outfiles = [self.make_trilegal_sfh(outfile=outfile_fmt % i,
                                           **mk_tri_sfh_kw)
                    for i in range(nsfhs)]

        return outfiles

    def compare_tri_match(self, trilegal_catalog, filter1, filter2,
                          outfig=None):
        '''
        Two plots, one M/H vs Age for match and trilegal, the other
        sfr for match vs age and number of stars of a given age for trilegal.
        '''
        sgal = rsp.Galaxies.simgalaxy(trilegal_catalog, filter1=filter1,
                                  filter2=filter2)
        sgal.lage = sgal.data.get_col('logAge')
        sgal.mh = sgal.data.get_col('[M/H]')
        issfr, = np.nonzero(self.sfr > 0)
        age_bins = np.digitize(sgal.lage, self.lagef[issfr])
        mean_mh= [np.mean(sgal.mh[age_bins==i]) for i in range(len(issfr))]

        bins = self.lagei
        sfr = np.array(np.histogram(sgal.lage, bins=bins)[0], dtype=float)

        fig, (ax1, ax2) = plt.subplots(figsize=(8,8), ncols=2, sharex=True)
        # should be density, weighted by number anyway..
        ax1.plot(sgal.lage, sgal.mh, '.', color='grey')
        ax1.plot(self.lagei[issfr], mean_mh, linestyle='steps', color='navy',
                 lw=3, label='TRILEGAL')

        ax1.plot(self.lagei[issfr], self.mh[issfr], linestyle='steps', lw=3,
                color='k', label='MATCH')
        ax1.fill_between(self.lagei[issfr],
                         self.mh[issfr] + self.mh_disp[issfr],
                         self.mh[issfr] - self.mh_disp[issfr],
                         lw=2, color='red', alpha=0.2)
        ax1.set_ylabel('$[M/H]$', fontsize=20)
        ax1.set_xlabel('$\log {\\rm Age (yr)}$', fontsize=20)
        ax1.legend(loc=0, frameon=False)

        ax2.plot(bins[:-1], sfr/(np.sum(sfr)), linestyle='steps', color='navy',
                lw=3, label='TRILEGAL')
        ax2.plot(self.lagei, self.sfr/np.sum(self.sfr),
                 linestyle='steps', lw=2, color='k', label='MATCH')
        ax2.set_ylabel('$ {\propto \\rm SFR}$', fontsize=20)
        ax2.set_xlabel('$\log {\\rm Age (yr)}$', fontsize=20)
        ax2.legend(loc=0, frameon=False)
        ax2.set_xlim(8, 10.5)
        if outfig is not None:
            fig.savefig(outfig, dpi=150)
