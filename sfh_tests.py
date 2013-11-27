import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import random
import logging

from astroML.stats import binned_statistic_2d
import brewer2mpl

import ResolvedStellarPops as rsp
import ResolvedStellarPops.convertz as convertz
from TPAGBparams import table_src, snap_src
import galaxy_tests


logging.basicConfig(filename='sfh_tests.log',level=logging.DEBUG)
logger = logging.getLogger()
logger.info('start of run')
angst_data = rsp.angst_tables.AngstTables()


def number_of_stars(gal=None, exclude_region='default', mag_below_trgb=2.,
                    galaxy_information=None, factor=2.):
    '''
    Count the number of rgb and agb stars (in F814W or F160W)

    ARGS:
    gal: rsp.Galaxies.galaxy instance
    galaxy_information: if no gal, dictonary with target, mag2,
                        filter2, filter1
    exclude_region: mags above and below the trgb to exclude from
                    star counts [default: factor times the trgb mag error]
    mag_below_trgb: mags below the trgb to call rgb stars [2]
    factor: factor times the trgb mag error to exclude

    RETURNS:
    nrgb: int number of rgb stars
    nagb: int number of agb stars
    exclude_dict: dict of trgb, trgb_err, and factor

    Use galaxy_information if its a trilegal cmd or galaxy instance
    for data.

    There are two reasons to use exclude_region.
    1) the TRGB mag is uncertain, and when counting agb stars,
       there should be no chance of scattering up rgb stars.
    2) TPAGB stars can be below the TRGB.

    TO DO:

    Adapt for recent star formation:
    Add color information or verts to the stars in region call.
    Will need a robust contamination code.
    '''

    # use with either data or simulations
    galaxy_information = galaxy_information or {}
    if gal is None:
        target = galaxy_information.get('target')
        mag2 = galaxy_information.get('mag2')
        filter2 = galaxy_information.get('filter2')
        filter1 = galaxy_information.get('filter1')
    else:
        target = gal.target
        mag2 = gal.mag2
        filter2 = gal.filter2
        filter1 = gal.filter1

    if filter2 == 'F814W':
        extra_key = '%s,%s' % (filter1, filter2)
        trgb = angst_data.get_item(target, 'mTRGB', extra_key)
        trgb_err = angst_data.get_item(target, 'mTRGB_err',
                                       extra_key=extra_key)
    if gal.filter2 == 'F160W':
        trgb = rsp.fileIO.item_from_row(angst_data.snap_tab3, 'target',
                                        target.upper(), 'mTRGB_F160W')
        trgb_err = rsp.fileIO.item_from_row(angst_data.snap_tab3,
                                            'target', target.upper(),
                                            'mTRGB_F160W_err')

    if exclude_region == 'default':
        exclude_region = trgb_err * factor
    else:
        # overwrite factor and use the exclude region as the mag above
        # and below. Then trgb_err must be set to 1 because it is always
        # multiplied by factor when used.
        trgb_err = 1.
        factor = exclude_region

    # rgb will go from mag_below_trgb to trgb + exclude_region
    nrgb = len(rsp.math_utils.between(mag2, trgb + mag_below_trgb,
                                      trgb + exclude_region))

    # agb will go from trgb - exclude_region to very bright
    nagb = len(rsp.math_utils.between(mag2, trgb - exclude_region, -99.))

    # add trgb_err to galaxy instance, could also add factor...
    if gal is not None:
        gal.trgb_err = trgb_err

    exclude_dict = {'trgb': trgb, 'trgb_err': trgb_err, 'factor': factor}

    return nrgb, nagb, exclude_dict


class AncientGalaxies(object):
    '''
    This class is to set up a run for VarySFHs.
    write_trgb_table: writes a table with trgb and the number of agb stars and
    rgb stars (using number_of_stars).
    number_of_stars: finds the number of rgb and agb stars in the data.
    read_trgb_table: reads the table write_trgb_table writes into recarray
    '''
    def __init__(self):
        pass

    def write_trgb_table(self, targets='ancients', mag_below_trgb=[2, 1.5],
                         exclude_region='default'):
        '''
        write the trgb, trgb_err, and the number of stars (opt, ir) rgb and agb
        '''
        if targets == 'ancients':
            if type(exclude_region) is list:
                ex_str = '_'.join((str(f) for f in exclude_region))
            else:
                ex_str = str(exclude_region)
            tstring = '%s_%s' % (targets, ex_str)
            targets = galaxy_tests.ancients()
        else:
            print 'need to write how to read other targets'

        fits_srcs = [snap_src + '/data/angst_no_trim', 'default']
        factors = []
        gal_dict = {}
        for i, band in enumerate(['opt', 'ir']):
            gals = rsp.Galaxies.galaxies([galaxy_tests.load_galaxy(t, band=band,
                                                      fits_src=fits_srcs[i])
                                          for t in targets])
            if type(exclude_region) is list:
                exreg = exclude_region[i]
            else:
                exreg = exclude_region

            for gal in gals.galaxies:
                if 'ngc404-deep' == gal.target:
                    gal.target = 'ngc404'
                else:
                    gal.target = gal.target.lower()

                nrgb, nagb, ex_dict = number_of_stars(gal,
                                                      mag_below_trgb=mag_below_trgb[i],
                                                      exclude_region=exreg)

                gdict = {'%s_trgb' % band: gal.trgb,
                         '%s_trgb_err' % band: ex_dict['trgb_err'],
                         'n%s_rgb' % band: nrgb,
                         'n%s_agb' % band: nagb}
                try:
                    gal_dict[gal.target].update(gdict)
                except KeyError:
                    gal_dict[gal.target] = gdict
            factors.append(ex_dict['factor'])
        fmt = '%(target)s '
        fmt += '%(opt_trgb).2f %(opt_trgb_err).2f %(nopt_rgb)i %(nopt_agb)i '
        fmt += '%(ir_trgb).2f %(ir_trgb_err).2f %(nir_rgb)i %(nir_agb)i \n'
        outfile = os.path.join(table_src, '%s_galaxies.dat' % tstring)

        header = '# mags below trgb: optical=%.2f nir=%.2f \n' % (mag_below_trgb[0],
                                                                  mag_below_trgb[1])
        if len(np.unique(factors)) > 1:
            header += '# excluded trgb region: %.2f * opt trgb mag_err \n' % factors[0]
            header += '# excluded trgb region: %.2f * ir trgb mag_err \n' % factors[1]
        else:
            header += '# excluded trgb region: %.2f * trgb mag_err \n' % factors[0]
        header += '# target opt_trgb opt_trgb_err nopt_rgb nopt_agb '
        header += 'nir_trgb nir_trgb_err nir_rgb nir_agb \n'

        with open(outfile, 'w') as out:
            out.write(header)
            for k, v in gal_dict.items():
                gal_dict[k]['target'] = k
                out.write(fmt % gal_dict[k])
        print 'wrote %s' % outfile
        return outfile

    def read_trgb_table(self, table_name):
        dtype = [('target', '|S16'), ('opt_trgb', '<f8'),
                 ('opt_trgb_err', '<f8'), ('nopt_rgb', '<f8'),
                 ('nopt_agb', '<f8'), ('ir_trgb', '<f8'),
                 ('ir_trgb_err', '<f8'), ('nir_rgb', '<f8'),
                 ('nir_agb', '<f8')]
        data = np.genfromtxt(table_name, dtype=dtype)
        self.data = data.view(np.recarray)

        # read the exclude factor around trgb_err
        with open(table_name, 'r') as f:
            lines = f.readlines()
        eline = [l for l in lines if 'excluded' in l]
        mline, = [l for l in lines if 'mags' in l]
        if len(eline) > 1:
            self.factor = np.array([el.split(':')[1].split('*')[0]
                                    for el in eline], dtype=float)
        else:
            factor = eline.split(':')[1].split('*')[0]
            self.factor = np.array([factor, factor], dtype=float)

        # read the mag_below offsets for each filter
        opt_off = float(mline.split('=')[1].split('nir')[0])
        nir_off = float(mline.strip().split('=')[-1])
        self.offsets = [opt_off, nir_off]


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
    if 'match' in file_origin.lower():
        data = rsp.match_utils.read_binned_sfh(filename)
    else:
        print 'please add a new data reader'
    
    if 'grid' in file_origin.lower():
        # at least have a uniform error that is 20% the smallest sfr.
        data.sfr_errp = np.min(data.sfr[data.sfr > 0]) * frac
        # now take 10% as nominal error in each sfr bin.
        data.sfr_errp[data.sfr > 0] = data.sfr[data.sfr > 0] * frac
        data.sfr_errm[data.sfr > 0] = data.sfr[data.sfr > 0] * frac
    return data


class StarFormationHistories(object):
    '''Make TRILEGAL star formation history files from MATCH'''
    def __init__(self, sfh_file, file_origin):
        self.base, self.name = os.path.split(sfh_file)
        self.data = parse_sfh_data(sfh_file, file_origin)

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
            print 'Warning: this method was designed for sfr, not [M/H]'
        # load in values this way in case I want to move this to its own
        # function
        val_arr = self.data.__getattribute__(attr)
        errm_arr = self.data.__getattribute__('%s_errm' % attr)
        errp_arr = self.data.__getattribute__('%s_errp' % attr)

        rand_arr = np.array([])
        # don't want negative sfr values. If not sfr, don't care.
        if attr == 'sfr':
            lowlim = 0
        else:
            lowlim = -np.inf

        for val, errm, errp in zip(val_arr, errm_arr, errp_arr):
            if errp == errm and errp > 0:
                # even uncertainties, easy.
                new_arr = np.random.normal(val, errp,  npoints)
            elif errp != 0 and errm != 0:
                # stitch two gaussians together
                pos_gauss = np.random.normal(val, errp,  npoints)
                neg_gauss = np.random.normal(val, errm, npoints)
                new_arr = np.concatenate([pos_gauss[pos_gauss >= val],
                                          neg_gauss[neg_gauss <= val]])

            elif errp == 0 and errm != 0:
                # no positive uncertainties
                neg_gauss = np.random.normal(val, errm, npoints)
                new_arr = neg_gauss[neg_gauss <= val]
            elif errp != 0 and errm == 0:
                # no negative uncertainties
                pos_gauss = np.random.normal(val, errp,  npoints)
                new_arr = pos_gauss[pos_gauss >= val]
            else:
                # um... no errors, why was this called
                print 'Warning: no uncertainties'
                new_arr = np.ones(4) * val
            new_arr = new_arr[new_arr > lowlim]
            rand_arr = np.append(rand_arr, random.choice(new_arr))
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

        if outfile == 'default':
          outfile = os.path.join(self.base,
                                 self.name.replace('.sfh', '.tri.dat'))
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
            # Choosing instead from within dispersion.
            interped_mh = self.interp_null_values()
            #mh = self.random_draw_within_uncertainty('mh')
            disp = np.median(self.data.mh_disp[np.nonzero(self.data.mh_disp)])/2.
            mh = np.array([random.choice(np.random.normal(interped_mh[i], disp,
                                                          2e5))
                           for i in range(len(interped_mh))])
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

    def plot_random_arrays(self, attr_str, val_arrs=None, ax=None, outfile=None,
                           from_files=False):
        '''
        val_arrs are random, attr_str is used to find the best fit values listed
        in the data table.
        after making a bunch of arrays that sample the sfr or mh uncertainties
        plot up where they are.
        '''
        if from_files is True:
            val_arrs = self.load_random_arrays(attr_str)

        assert val_arrs is not None, 'either specify val_arrs or set from_files'

        val_arr = self.data.__getattribute__(attr_str)
        errm_arr = self.data.__getattribute__('%s_errm' % attr_str)
        errp_arr = self.data.__getattribute__('%s_errp' % attr_str)

        if 'sfr' in attr_str:
            ylab = '${\\rm SFR}$'
        elif 'm' in attr_str:
            ylab = '${\\rm [M/H]}$'
        elif 'fe' in attr_str:
            ylab = '${\\rm [Fe/H]}$'

        if ax is None:
            fig, ax = plt.subplots(figsize=(12,12))

        [ax.errorbar(self.data.lagei, val_arrs[i], linestyle='steps-mid',
                     color='k', alpha=0.3) for i in range(len(val_arrs))]
        ax.errorbar(self.data.lagei, val_arr, [errm_arr, errp_arr],
                    linestyle='steps-mid', lw=2, color='r')
        ax.set_ylabel(ylab, fontsize=20)
        ax.set_xlabel('$\log {\\rm Age (yr)}$', fontsize=20)
        ax.set_xlim(8, 10.5)
        if outfile is not None:
            plt.savefig(outfile, dpi=300)
        return ax

    def make_many_trilegal_sfhs(self, nsfhs=100, mk_tri_sfh_kw=None):
        '''
        make nsfhs number of trilegal sfh input files.
        '''
        mk_tri_sfh_kw = mk_tri_sfh_kw or {}

        if nsfhs > 1:
            self.mc = True

        outfile_fmt = mk_tri_sfh_kw.get('outfile_fmt', None)
        
        if outfile_fmt is None:
            new_dir = os.path.join(self.base, 'mc/')
            rsp.fileIO.ensure_dir(new_dir)
            outfile_fmt = os.path.join(new_dir,
                                       os.path.split(self.name)[1].replace('.zc.sfh',
                                                                           '_%02i.tri.dat'))
        else:
            # need to iterate outside of dict... see below.
            del mk_tri_sfh_kw['outfile_fmt']

        # update any passed kw to make_trilegal_sfh with defaults.
        mk_tri_sfh_kw = dict({'random_sfr': True, 'random_z': True}.items()
                             + mk_tri_sfh_kw.items())
        
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
        # should be density, weighted by number anyway...
        ax1.plot(sgal.lage, sgal.mh, '.', color='grey')
        ax1.plot(self.lagei[issfr], mean_mh, linestyle='steps', color='navy',
                 lw=3, label='TRILEGAL')

        ax1.plot(self.lagei[issfr], self.mh[issfr], linestyle='steps', lw=3,
                color='k', label='MATCH')
        ax1.fill_between(self.lagei[issfr], self.mh[issfr] + self.mh_disp[issfr],
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
            fig.savefig(outfig, dpi=300)


class VarySFHs(StarFormationHistories, AncientGalaxies):
    '''
    run several variations of the age sfr z from MATCH SFH to produce
    simulated CMDs, LFs, and nagb/nrgb ratios.

    can make plots, and save summary files.

    First, use AncientGalaxies to create a table of trgb and nrgb and
    nagb stars from data and pass this as the table_file.

    NOTE: AncientGalaxies and just takes all the stars within some mag
    below trgb. Something better is needed for recent SF.

    The main method is vary_the_SFH, which sets up the files, calls
    StarFormationHistories and calls the other methods to write out the
    LFs and ratios.
    '''
    def __init__(self, galaxy_input, sfh_file, file_origin, target=None,
                 table_file='default', ast=False):
        StarFormationHistories.__init__(self, sfh_file, file_origin)
        self.galaxy_input = galaxy_input

        if target is None:
            gname = os.path.split(self.galaxy_input)[1]
            target = gname.split('_')[1].replace('.dat', '').lower()
        self.target = target

        if table_file is 'default':
            table_file = os.path.join(table_src, 'ancients_galaxies.dat')
        else:
            print 'reading from table %s', table_file
        # Reads in the data as well as the mag offsets and factor in the
        # exclude region.
        self.ags = AncientGalaxies()
        self.ags.read_trgb_table(table_file)

        # load stars in data RGB and AGB region as well as TRGBs.
        # this is just a call to self.ags.data
        self.load_data_for_normalization()

        # should we use ast corrections?
        self.ast = ast
        if self.ast is False:
            print 'not using artificial star corrections'

    def prepare_trilegal_sfr(self, make_many_kw=None):
        '''call make_many_trilegal_sfhs'''
        make_many_kw = make_many_kw or {}
        self.sfr_files = self.make_many_trilegal_sfhs(**make_many_kw)

    def prepare_galaxy_input(self, object_mass=None, dry_run=False):
        '''
        write the galaxy input file from a previously written template.
        simply overwrites the filename line to link to the new sfr
        file.
        '''
        self.galaxy_inputs = []

        lines = open(self.galaxy_input).readlines()
        # line that links to sfr file.
        extra = ' '.join(lines[-3].split(' ')[1:])

        if object_mass is not None:
            extra2 = ' '.join(lines[-6].split()[1:])

        for i, sfr_file in enumerate(self.sfr_files):
            lines[-3] = ' '.join([self.sfr_files[i], extra])
            if object_mass is not None:
                lines[-6] = ' '.join(['%.4e' % object_mass, extra2])
            new_out = os.path.join(self.outfile_loc,
                                   os.path.split(self.galaxy_input)[1].replace('.dat', '_%003i.dat' % i))
            if dry_run is False:
                with open(new_out, 'w') as f:
                    f.write(''.join(lines))
                print 'wrote %s' % new_out
            self.galaxy_inputs.append(new_out)

    def prepare_outfiles(self, outfile_loc='default', extra_directory=None,
                         clean_first=False):
        '''
        prepare outfiles for vary_the_sfh.
        opt_lf, ir_lf
        opt_lf_noagb, ir_lf_noagb
        opt_mass_met, ir_mass_met
        narratio

        NOTE: all attrs saved are file objects.
        '''
        if outfile_loc == 'default':
            if extra_directory is None:
                self.outfile_loc = os.path.join(snap_src, 'models', 'varysfh',
                                                self.target, 'mc/')
                extra_directory = ''
            else:
                self.outfile_loc = os.path.join(snap_src, 'models', 'varysfh',
                                                self.target, extra_directory,
                                                'mc/')
                extra_directory += '_'

        rsp.fileIO.ensure_dir(self.outfile_loc)
        if clean_first is True:
            self.remove_files()

        # N agb/rgb ratio file
        hfmt =  '%(target)s %(nopt_rgb)i %(nopt_agb)i %(nir_rgb)i '
        hfmt += '%(nir_agb)i %(opt_ar_ratio).3f %(ir_ar_ratio).3f '
        hfmt += '%(opt_ar_ratio_err).3f  %(ir_ar_ratio_err).3f \n'

        header = '# target nopt_rgb nopt_agb nir_rgb nir_agb opt_ar_ratio '
        header += 'ir_ar_ratio opt_ar_ratio_err ir_ar_ratio_err \n'

        fnames = ['opt_lf', 'ir_lf', 'opt_lf_noagb', 'ir_lf_noagb', 'narratio',
                  'opt_mass_met', 'ir_mass_met']

        for fname in fnames:
            name_fmt = '%s%s_%s.dat' % (extra_directory, self.target, fname)
            name =  os.path.join(self.outfile_loc, name_fmt)
            exists = os.path.isfile(name)
            self.__setattr__('%s_file' % fname, open(name, 'a'))
            # it's ugly to keep writing the header.
            self.narratio_fmt = hfmt
            if exists is False and fname == 'narratio':
                self.narratio_file.write(header)
        
        # how to name the trilegal sfr files
        new_fmt = self.name + '_tri_%003i.sfr'
        self.sfr_outfilefmt = os.path.join(self.outfile_loc, new_fmt)
        return
    
    def close_files(self):
        '''close all open files'''
        [f.close() for f in self.__dict__.values() if type(f) is file]

    def write_LF(self, opt_norm, ir_norm, opt_file, ir_file, hist_it_up=True,
                 inds=None):
        '''
        write the LF to files.
        ARGS:
        opt_file, ir_file: write out files, can be file objects.
        opt_rec and ir_rec are indices of sgal f814w f160w.
        hist_it_up: use galaxy_tests.hist_it_up or bayesian_blocks
        NOTE: must have attr sgal loaded.
        Formatting is a repeat of first line is histogram, next line is bins[1:]
        '''
        inds = inds or np.arange(len(self.sgal.mag2))

        # load mags, set here in case incase I want to pass this in the future.
        opt_mag, ir_mag = self.load_trilegal_data()

        if hist_it_up is True:
            # hist it up!
            opt_hist, opt_bins = galaxy_tests.hist_it_up(opt_mag, threash=5)
            ir_hist, ir_bins = galaxy_tests.hist_it_up(ir_mag, threash=5)
        else:
            nbins = (np.max(opt_mag) - np.min(opt_mag)) / 0.1
            opt_hist, opt_bins = np.histogram(opt_mag, bins=nbins)

            nbins = (np.max(ir_mag) - np.min(ir_mag)) / 0.1
            ir_hist, ir_bins = np.histogram(ir_mag, bins=nbins)

        # scale the simulated LF to match the data LF
        opt_hist = np.array(opt_hist, dtype=float) * opt_norm
        ir_hist = np.array(ir_hist, dtype=float) * ir_norm

        np.savetxt(opt_file, (opt_hist, opt_bins[1:]), fmt='%g')
        np.savetxt(ir_file, (ir_hist, ir_bins[1:]), fmt='%g')

        print 'wrote %s, %s' % (opt_file.name, ir_file.name)
        return

    def count_stars_from_hist(self, opt_hist, opt_bins, ir_hist, ir_bins):
        ratio_data = {}
        for i, (hist, bins, band) in enumerate(zip([opt_hist, ir_hist],
                                                 [opt_bins, ir_bins],
                                                 ['opt', 'ir'])):
            trgb = self.__getattribute__('%s_trgb' % band)
            trgb_err = self.__getattribute__('%s_trgb_err' % band)
            norm = self.ags.offsets[i]
            irgb = rsp.math_utils.between(bins, trgb + norm,
                                          trgb + trgb_err * self.ags.factor[i])
            iagb = rsp.math_utils.between(bins,
                                          trgb - trgb_err  * self.ags.factor[i], 10.)

            nrgb = np.sum(hist[irgb])
            nagb = np.sum(hist[iagb])
            ratio_data['%s_ar_ratio' % band] = nagb / nrgb
            ratio_data['%s_ar_ratio_err' % band] = galaxy_tests.count_uncert_ratio(nagb, nrgb)
            ratio_data['n%s_rgb' % band] = nrgb
            ratio_data['n%s_agb'% band] = nagb
        return ratio_data

    def write_ratio(self, opt_rgb, opt_agb, ir_rgb, ir_agb, out_file):
        '''
        write the numbers of stars, the ratio, and the poisson
        uncertainty in the ratio.
        see self.prepare_outfiles for file format information

        ARGS:
        opt_rgb: list of optical rgb star indices
        opt_agb: list of optical agb star indices
        ir_rgb: list of NIR rgb star indices
        ir_agb: list of NIR agb star indices
        out_file: is an opened file object.
        '''

        nopt_rgb = float(len(opt_rgb))
        nopt_agb = float(len(opt_agb))
        nir_rgb = float(len(ir_rgb))
        nir_agb = float(len(ir_agb))
        out_dict = {'target': self.target,
                    'opt_ar_ratio': nopt_agb / nopt_rgb,
                    'ir_ar_ratio': nir_agb / nir_rgb,
                    'opt_ar_ratio_err': galaxy_tests.count_uncert_ratio(nopt_agb, nopt_rgb),
                    'ir_ar_ratio_err': galaxy_tests.count_uncert_ratio(nir_agb, nir_rgb),
                    'nopt_rgb': nopt_rgb,
                    'nopt_agb': nopt_agb,
                    'nir_rgb': nir_rgb,
                    'nir_agb': nir_agb}

        out_file.write(self.narratio_fmt % out_dict)
        print 'wrote %s' % out_file.name

    def write_mass_met_file(self):
        '''
        Write the magnitudes, masses, and [M/H] values of tpagb stars
        to files.

        Format of the file is blocks of three lines, mag, mass, [M/H] for
        each sfh run.
        '''
        # load tpagb indices
        self.sgal.all_stages('TPAGB')
        if self.ast is True:
            # recovered tpagb stars
            opt_inds = np.intersect1d(self.sgal.itpagb, self.sgal.rec)
            ir_inds = np.intersect1d(self.sgal.itpagb, self.sgal.ir_rec)
        else:
            opt_inds = self.sgal.itpagb
            ir_inds = self.sgal.itpagb

        mag2 = self.sgal.mag2
        mag4 = self.sgal.data.get_col('F160W')
        for mag, inds, ofile in zip([mag2, mag4], [opt_inds, ir_inds],
                                    [self.opt_mass_met_file,
                                     self.ir_mass_met_file]):
            mass = self.sgal.data.get_col('m_ini')[inds]
            mag = mag[inds]
            mh = self.sgal.data.get_col('[M/H]')[inds]
            np.savetxt(ofile, (mag, mass, mh), fmt='%g')
            
        print 'wrote %s, %s' % (self.opt_mass_met_file.name,
                                self.ir_mass_met_file.name)

    def load_data_for_normalization(self):
        '''load the numbers of data rgb and agb stars from self.ags'''
        self.nopt_rgb = rsp.fileIO.item_from_row(self.ags.data, 'target',
                                                 self.target, 'nopt_rgb')
        self.nir_rgb = rsp.fileIO.item_from_row(self.ags.data, 'target',
                                                self.target, 'nir_rgb')
        self.opt_trgb = rsp.fileIO.item_from_row(self.ags.data, 'target',
                                                self.target, 'opt_trgb')
        self.ir_trgb = rsp.fileIO.item_from_row(self.ags.data, 'target',
                                                self.target, 'ir_trgb')
        self.opt_trgb_err = rsp.fileIO.item_from_row(self.ags.data, 'target',
                                                self.target, 'opt_trgb_err')
        self.ir_trgb_err = rsp.fileIO.item_from_row(self.ags.data, 'target',
                                                self.target, 'ir_trgb_err')
        return self.nopt_rgb, self.nir_rgb

    def load_asts(self):
        '''load rsp.Galaxies.artificial_star_tests objects'''
        fake_files = galaxy_tests.get_fake_files(self.target.upper())
        if self.ast is True:
            ast_objs = [rsp.Galaxies.artificial_star_tests(f)
                        for f in fake_files]
            self.ast_objs = ast_objs

    def read_trilegal_catalog(self, trilegal_output, filter1='F606W'):
        '''read the trilegal cat and does ast corrections.'''
        sgal = rsp.Galaxies.simgalaxy(trilegal_output, filter1=filter1,
                                      filter2='F814W')
        if self.ast is True and not hasattr(self, 'ast_objs'):
            # should be loaded outside this method if mc is running.
            self.load_asts()

            # ast correction
            rsp.Galaxies.ast_correct_trilegal_sim(sgal, asts_obj=self.ast_objs)
            sgal.load_ast_corrections()
        
        self.sgal = sgal

    def shift_ir(self):
        '''shift ir mags to they agree with opt trgb'''

        trgb_color = angst_data.get_item(self.target, 'mean_color')

        opt_mag = self.sgal.data.get_col('F814W')
        opt_mag1 = self.sgal.data.get_col(self.filter1)
        opt_color = opt_mag1 - opt_mag
        
        ir_mag = self.sgal.data.get_col('F160W')
        # closest star in the trilegal catalog to star on the trgb
        ind, dist = rsp.math_utils.min_dist2d(trgb_color, self.opt_trgb,
                                              opt_color, opt_mag)
        # correction for ir mag
        offset = self.ir_trgb - ir_mag[ind]
        return offset

    def load_trilegal_data(self):
        '''load trilegal F814W and F160W mags'''
        extra = ''
        if self.ast is True:
            extra = '_cor'
        print extra, 'extra'
        opt_mag = self.sgal.data.get_col('%s%s' % ('F814W', extra))
        ir_mag = self.sgal.data.get_col('%s%s' % ('F160W', extra))
        offset = self.shift_ir()
        print 'IR OFFSET:', offset
        ir_mag += offset
        # it's a matter of memory to not keep all the non-recovered stars
        # the hard coded limits just help not having to set ylims in LF plots.
        if self.ast is True:
            opt_mag = opt_mag[np.isfinite(opt_mag)]
            opt_mag = opt_mag[opt_mag < 29]
            ir_mag = ir_mag[np.isfinite(ir_mag)]
            ir_mag = ir_mag[ir_mag < 26]
        return opt_mag, ir_mag

    def do_normalization(self, filter1=None, trilegal_output=None,
                         hist_it_up=True):
        '''Do the normalization and save small part of outputs.'''
        if not hasattr(self, 'sgal'):
            assert trilegal_output is not None, \
                'need sgal loaded or pass trilegal catalog file name'
        if self.mc is True:
            self.read_trilegal_catalog(trilegal_output, filter1=filter1)

        opt_mag, ir_mag = self.load_trilegal_data()

        # Recovered stars in simulated RGB region.
        sopt_rgb = self.sgal.stars_in_region(opt_mag,
                                             self.opt_trgb + self.ags.offsets[0],
                                             self.opt_trgb + \
                                             self.opt_trgb_err * \
                                             self.ags.factor[0])
        sir_rgb = self.sgal.stars_in_region(ir_mag,
                                            self.ir_trgb + self.ags.offsets[1],
                                            self.ir_trgb + \
                                            self.ir_trgb_err * self.ags.factor[1])

        # Recovered stars in simulated AGB region.
        sopt_agb = self.sgal.stars_in_region(opt_mag, self.opt_trgb - \
                                             self.opt_trgb_err * \
                                             self.ags.factor[0], 10.)
        sir_agb = self.sgal.stars_in_region(ir_mag, self.ir_trgb - \
                                            self.ir_trgb_err  * \
                                            self.ags.factor[1], 10.)

        # normalization
        opt_norm = self.nopt_rgb / float(len(sopt_rgb))
        ir_norm = self.nir_rgb / float(len(sir_rgb))

        # for tracking non agb stars
        itpagb = self.sgal.stage_inds('TPAGB')
        non_tpagb = list(set(np.arange(len(self.sgal.mag2))) - set(itpagb))

        # random sample the data distribution
        rands = np.random.random(len(opt_mag))
        opt_ind, = np.nonzero(rands < opt_norm)
        rands = np.random.random(len(ir_mag))
        ir_ind, = np.nonzero(rands < ir_norm)

        # scaled rgb: ast + norm + in rgb
        opt_rgb = list(set(opt_ind) & set(sopt_rgb))
        ir_rgb = list(set(ir_ind) & set(sir_rgb))

        #print len(opt_rgb), self.nopt_rgb
        #print len(ir_rgb), self.nir_rgb

        # scaled agb
        opt_agb = list(set(opt_ind) & set(sopt_agb))
        ir_agb = list(set(ir_ind) & set(sir_agb))

        self.opt_norm = opt_norm
        self.ir_norm = ir_norm

        #save LF in both filters
        #if self.mc is True:
        # only need to save the LF and the ratios.
        self.write_mass_met_file()
        self.write_LF(opt_norm, ir_norm, self.opt_lf_file, self.ir_lf_file,
                      hist_it_up=hist_it_up)
        self.write_LF(opt_norm, ir_norm, self.opt_lf_noagb_file, self.ir_lf_noagb_file,
                      hist_it_up=hist_it_up, inds=non_tpagb)
        self.write_ratio(opt_rgb, opt_agb, ir_rgb, ir_agb,
                         self.narratio_file)

        return opt_norm, ir_norm, opt_ind, ir_ind

    def vary_the_SFH(self, cmd_input_file, make_many_kw=None, dry_run=False,
                     diag_plots=False, hist_it_up=True,
                     extra_directory='default', clean_first=False,
                     add_stage_lfs=None, filter1=None):
        '''
        make the sfhs, make the galaxy inputs, run trilegal. For no trilegal
        runs, set dry_run True.
        '''
        agb_model = cmd_input_file.replace('cmd_input_', '').lower()
        if extra_directory == 'default':
            agb_mod = agb_model.split('.')[0]
        else:
            agb_mod = None
        self.agb_mod = agb_mod

        self.mc = False
        nsfhs = make_many_kw.get('nsfhs', 50)
        if nsfhs > 1:
            self.mc = True

        self.prepare_outfiles(extra_directory=agb_mod, clean_first=clean_first)

        make_many_kw = make_many_kw or {}
        if not 'mk_tri_sfh_kw' in make_many_kw.keys():
            make_many_kw['mk_tri_sfh_kw'] = {}

        make_many_kw['mk_tri_sfh_kw'] = dict({'outfile_fmt':
            self.sfr_outfilefmt, 'dry_run': dry_run}.items() +
            make_many_kw['mk_tri_sfh_kw'].items())

        self.prepare_trilegal_sfr(make_many_kw=make_many_kw)

        self.prepare_galaxy_input(dry_run=dry_run)

        tname = 'output_%s_%s' % (self.target, agb_model)
        trilegal_output = os.path.join(self.outfile_loc, tname)

        if self.ast is True:
            self.load_asts()
            filter1 = self.ast_objs[0].filter1
        else:
            filter1 = get_filter1(self.target)
            self.filter1 = filter1

        for galaxy_input in self.galaxy_inputs:
            rsp.TrilegalUtils.run_trilegal(cmd_input_file, galaxy_input,
                                           trilegal_output, rmfiles=False,
                                           dry_run=dry_run)
            self.do_normalization(trilegal_output=trilegal_output,
                                  filter1=filter1, hist_it_up=hist_it_up)
        self.close_files()
        if diag_plots is True:
            outfile_fmt = os.path.join(self.outfile_loc, '%s_random' % self.target)
            outfile_fmt += '_%s.png'
            [self.plot_random_arrays(attr_str, from_files=True,
                                     outfile=outfile_fmt % attr_str)
             for attr_str in ['sfr', 'mh']]
            
            self.compare_to_gal(hist_it_up=hist_it_up, add_stage_lfs=add_stage_lfs)

            self.plot_mass_met_table(self.opt_mass_met_file.name,
                                     self.ir_mass_met_file.name)
    
    def plot_lf_file(self, opt_lf_file, ir_lf_file, axs=None, plt_kw=None):
        '''needs work, but: plot the lf files.'''
        plt_kw = plt_kw or {}
        plt_kw = dict({'linestyle': 'steps',
                       'color': 'black',
                       'alpha': 0.2}.items() + plt_kw.items())
        if axs is None:
            fig, (axs) = plt.subplots(ncols=2, figsize=(12, 6))
            plt.subplots_adjust(right=0.95, left=0.05, wspace=0.1)

        for i, lf_file in enumerate([opt_lf_file, ir_lf_file]):

            with open(lf_file, 'r') as lff:
                lines = [l.strip() for l in lff.readlines()
                         if not l.startswith('#')]

            hists = [np.array(l.split(), dtype=float) for l in lines[0::2]]
            binss = [np.array(l.split(), dtype=float) for l in lines[1::2]]

            for hist, bins in zip(hists, binss):
                if len(hist) != len(bins):
                    continue
                axs[i].plot(bins, hist, **plt_kw)
        return axs

    def compare_to_gal(self, hist_it_up=True, narratio=True, no_agb=False,
                       add_stage_lfs=None):
        '''
        Plot the LFs and galaxy LF.

        ARGS:
        hist_it_up: Use hist_it_up or bayseyn blocks
        narratio: overlay NRGB, NAGB, and NAGB/NRGB +/- err
        no_agb: plot the LF without AGB stars
        add_stage_lfs: (list) add LF of specific stages

        RETURNS:
        ax1, ax2: axes instances created for the plot.

        '''
        # load galaxy data
        ir_gal = galaxy_tests.load_galaxy(self.target, band='ir')
        fits_src = snap_src + '/data/angst_no_trim'
        opt_gal = galaxy_tests.load_galaxy(self.target, band='opt',
                                           fits_src=fits_src)
        # make galaxy histograms
        if hist_it_up is True:
            opt_hist, opt_bins = galaxy_tests.hist_it_up(opt_gal.mag2,
                                                        threash=5)
            ir_hist, ir_bins = galaxy_tests.hist_it_up(ir_gal.mag2, threash=5)
        else:
            nbins = (np.max(opt_gal.mag2) - np.min(opt_gal.mag2)) / 0.1
            opt_hist, opt_bins = np.histogram(opt_gal.mag2, bins=nbins)

            nbins = (np.max(ir_gal.mag2) - np.min(ir_gal.mag2)) / 0.1
            ir_hist, ir_bins = np.histogram(ir_gal.mag2, bins=nbins)

        if self.mc is True:
            # plot lfs from simulations (and initialize figure)
            ax1, ax2 = self.plot_lf_file(self.opt_lf_file.name,
                                         self.ir_lf_file.name)
        else:
            fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(12, 6), sharey=False)
            plt.subplots_adjust(right=0.95, left=0.1, wspace=0.1)

        # plot lfs without tpagb stars
        if no_agb is True:
            ax1, ax2 = self.plot_lf_file(self.opt_lf_noagb_file.name,
                                         self.ir_lf_noagb_file.name,
                                         axs=(ax1, ax2),
                                         plt_kw={'color': 'darkred', 'lw': 2,
                                                 'alpha': 1})
        if add_stage_lfs is not None and self.mc is False:
            if add_stage_lfs == 'all':
                add_stage_lfs = ['PMS', 'MS', 'SUBGIANT', 'RGB', 'HEB', 'RHEB',
                                 'BHEB', 'EAGB', 'TPAGB', 'POSTAGB', 'WD']
            if add_stage_lfs == 'default':
                add_stage_lfs = ['RGB', 'EAGB', 'TPAGB']

            nstages = len(add_stage_lfs)
            stage_lf_kw = {'linestyle': 'steps', 'lw': 2}
            if nstages < 3:
                cmap = brewer2mpl.get_map('Paired', 'Qualitative', 3)
                cols = cmap.mpl_colors[0::2]
            else:
                cmap = brewer2mpl.get_map('Paired', 'Qualitative', nstages)
                cols = cmap.mpl_colors

            sopt_mag, sir_mag = self.load_trilegal_data()
            for ax, mag, norm in zip([ax1, ax2], [sopt_mag, sir_mag],
                                     [self.opt_norm, self.ir_norm]):
                self.sgal.make_lf(mag, stages=add_stage_lfs,
                                  hist_it_up=hist_it_up)
                for i in range(nstages):
                    try:
                        hist = self.sgal.__getattribute__('i%s_lfhist' % add_stage_lfs[i].lower())
                    except AttributeError:
                        continue
                    bins = self.sgal.__getattribute__('i%s_lfbins' % add_stage_lfs[i].lower())
                    stage_lf_kw['color'] = cols[i]
                    stage_lf_kw['label'] = '$%s$' % add_stage_lfs[i]
                    ax.plot(bins[:-1], hist*norm, **stage_lf_kw)

            sopt_hist, sopt_bins = self.sgal.make_lf(sopt_mag,
                                                     hist_it_up=hist_it_up)
            sir_hist, sir_bins = self.sgal.make_lf(sir_mag,
                                                   hist_it_up=hist_it_up)

            sopt_hist = sopt_hist * self.opt_norm
            sir_hist = sir_hist * self.ir_norm
            stage_lf_kw['color'] = 'grey'
            lab = self.agb_mod or ''
            if lab != '':
                lab = '$%s$' % '\ '.join(lab.split('_')[1:])
            stage_lf_kw['label'] = lab
            ax1.plot(sopt_bins[:-1], sopt_hist, **stage_lf_kw)
            ax2.plot(sir_bins[:-1], sir_hist, **stage_lf_kw)
        # plot galaxy data
        dplot_kw = {'linestyle': 'steps', 'color': 'black', 'lw': 2,
                    'label': '$%s$' % self.target.upper()}
        ax2.plot(ir_bins[:-1], ir_hist, **dplot_kw)
        ax1.plot(opt_bins[:-1], opt_hist, **dplot_kw)

        # initialize add numbers to the plot
        if narratio is True:
            if add_stage_lfs is None or self.mc is True:
                # count stars from the saved file (which counts random asts)
                ratio_data = rsp.fileIO.readfile(self.narratio_file.name)
            else:
                # count stars from the loaded histograms
                ratio_data = self.count_stars_from_hist(sopt_hist, sopt_bins,
                                                        sir_hist, sir_bins)
            stext_kw = {'color': 'darkred', 'fontsize': 16, 'ha': 'center'}
            dtext_kw = {'color': dplot_kw['color'], 'fontsize': 16,
                        'ha': 'center'}

        for i, (ax, gal, offset, trgb_err, band) in enumerate(zip([ax1, ax2],
                                                   [opt_gal, ir_gal],
                                                   self.ags.offsets,
                                                   [self.opt_trgb_err,
                                                   self.ir_trgb_err],
                                                   ['opt', 'ir'])):
            ax.set_yscale('log')
            ax.set_ylim(1, ax.get_ylim()[1])
            ax.set_xlim(np.round(self.sgal.data.get_col('%s' % gal.filter2).min()), gal.comp50mag2)

            # vertical lines around the trgb exclude region
            ax.vlines(gal.trgb + trgb_err * self.ags.factor[i], *ax.get_ylim(),
                      lw=2, color='darkred')
            ax.vlines(gal.trgb - trgb_err * self.ags.factor[i], *ax.get_ylim(),
                      lw=2, color='darkred')
            # trgb line
            ax.vlines(gal.trgb, *ax.get_ylim(),
                      lw=2, color='black')
            # line at dim limit for rgb normalization
            ax.vlines(gal.trgb + offset, *ax.get_ylim(), lw=2, color='darkred')

            ax.legend(loc=1, frameon=False)
            ax.set_xlabel('$%s$' % gal.filter2, fontsize=20)

            if narratio is True:
                yval = 1.2  # text yloc found by eye, depends on fontsize

                # text xlocs are chosen in the middle of each region
                xrgb_val = np.mean([gal.trgb, gal.trgb + offset])
                xagb_val = np.mean([ax.get_xlim()[0], gal.trgb])
                xratio_val = np.mean([gal.trgb + offset, ax.get_xlim()[1]])

                # simulated nrgb and nagb are the mean values
                srgb_text = '$%i$' % np.mean(ratio_data['n%s_rgb' % band])
                sagb_text = '$%i$' % np.mean(ratio_data['n%s_agb' % band])

                # one could argue taking the mean isn't the best idea for
                # the ratio errors.
                sratio_text = '$%.3f\pm%.3f$' % (np.mean(ratio_data['%s_ar_ratio' % band]),
                                                 np.mean(ratio_data['%s_ar_ratio_err' % band]))

                # need to load the data nrgb and nagb, calculate the ratio
                # and error.
                nrgb = rsp.fileIO.item_from_row(self.ags.data, 'target',
                                                self.target, 'n%s_rgb' % band)
                nagb = rsp.fileIO.item_from_row(self.ags.data, 'target',
                                                self.target, 'n%s_agb' % band)
                dratio = nagb / nrgb
                dratio_err = galaxy_tests.count_uncert_ratio(nagb, nrgb)

                drgb_text = '$%i$' % nrgb
                dagb_text = '$%i$' % nagb
                dratio_text =  '$%.3f\pm%.3f$' % (dratio, dratio_err)

                ax.text(xrgb_val, yval, srgb_text, **stext_kw)
                ax.text(xagb_val, yval, sagb_text, **stext_kw)
                ax.text(xratio_val, yval, sratio_text, **stext_kw)

                yval += .7  # stack the text
                ax.text(xrgb_val, yval, drgb_text, **dtext_kw)
                ax.text(xagb_val, yval, dagb_text,  **dtext_kw)
                ax.text(xratio_val, yval, dratio_text, **dtext_kw)

        ax1.set_ylabel('$\#$', fontsize=20)
        plt.tick_params(labelsize=16)
        outfile = os.path.join(self.outfile_loc, '%s_lfs.png' % self.target)
        plt.savefig(outfile, dpi=300)
        print 'wrote %s' % outfile
        return ax1, ax2

    def plot_mass_met_table(self, opt_mass_met_file, ir_mass_met_file):
        fig = plt.figure(figsize=(9, 9))
        grid = ImageGrid(fig, 111,
                         nrows_ncols=(2, 2),
                         axes_pad=.5,
                         add_all=True,
                         label_mode="all",
                         cbar_location="top",
                         cbar_mode="each",
                         cbar_size="7%",
                         cbar_pad="2%",
                         aspect=0)
        cmaps = [plt.cm.get_cmap('jet', 9), plt.cm.gray_r]
        #cmap = 
        #cmap.set_bad('w', 1.)
        #fig, (axs) = plt.subplots(ncols=2, figsize=(8, 8), sharey=True)
        types = ['mean', 'count']
        k =-1
        for j in range(len(types)):
            for i, mass_met in enumerate([opt_mass_met_file,
                                          ir_mass_met_file]):
                k += 1            
                with open(mass_met, 'r') as mmf:
                    lines = [l.strip() for l in mmf.readlines()
                             if not l.startswith('#')]
    
                mag = np.concatenate([np.array(l.split(), dtype=float)
                                      for l in lines[0::3]])
                mass = np.concatenate([np.array(l.split(), dtype=float)
                                       for l in lines[1::3]])
                mh = np.concatenate([np.array(l.split(), dtype=float)
                                     for l in lines[2::3]])
                
                N, xedges, yedges = binned_statistic_2d(mag, mass, mh,
                                                        types[j], bins=50)
                im = grid[k].imshow(N.T, origin='lower',
                               extent=[xedges[0], xedges[-1], yedges[0],
                                       yedges[-1]],
                               aspect='auto', interpolation='nearest',
                               cmap=cmaps[j])
                grid[k].cax.colorbar(im)
                #grid[i].cax.set_label('$[M/H]$')

        grid.axes_all[0].set_ylabel('${\\rm Mass}\ (M_\odot)$', fontsize=20)
        grid.axes_all[2].set_ylabel('${\\rm Mass}\ (M_\odot)$', fontsize=20)
        grid.axes_all[2].set_xlabel('$F814W$', fontsize=20)
        grid.axes_all[3].set_xlabel('$F160W$', fontsize=20)
        fig.suptitle('$%s$' % self.target, fontsize=20)
        plt.savefig(os.path.join(self.outfile_loc,
                                '%s_mass_met.png' % self.target), dpi=300)
        return grid

    def remove_files(self):
        print 'backing up everything in %s first' % self.outfile_loc
        files = rsp.fileIO.get_files(self.outfile_loc, '*.*')
        bkdir = os.path.join(self.outfile_loc, 'bkup/')
        if os.path.isdir(bkdir):
            old_files = os.listdir(bkdir)
            if len(old_files) > 0:
                print 'overwriting a previous back up.'
        else:
            rsp.fileIO.ensure_dir(bkdir)
        [os.system('mv %s %s' % (f, bkdir)) for f in files]


def get_filter1(target, fits_src=None):
    '''get optical filter1'''
    fits_src = fits_src or snap_src + '/data/angst_no_trim'
    fitsfile, = rsp.fileIO.get_files(fits_src, '*%s*.fits' % target.lower())
    filter1 = fitsfile.split(target)[1].split('_')[1].split('-')[0].upper()
    return filter1


def convert_match_to_trilegal_sfh(sfr_dir, fileorigin='match',
                                  search_str='*.sfh',
                                  make_trilegal_sfh_kw=None):
    '''
    call StarFormationHistories.make_trilegal_sfh for a number of files in a
    directory.

    ARGS:
    sfr_dir     base location of match sfr files
    fileorigin  'match' sfr fileorigin, match or match-grid
    search_str  '*.sfh' extenstion of match sfh file
    make_trilegal_sfh_kw    {'random_sfr': False, 'random_z': False} kwargs to
                            send to make_trilegal_sfh
    RETURNS:
    list of trilegal SFR files written
    '''
    # make_trilegal_sfh kwargs will set both randoms to False by default
    make_trilegal_sfh_kw = make_trilegal_sfh_kw or {}
    make_trilegal_sfh_kw = dict({'random_sfr': False, 'random_z': False}.items() +
                                 make_trilegal_sfh_kw.items())
    # find match sfh files`
    sfhs = rsp.fileIO.get_files(sfr_dir, search_str)

    # load SFHs
    SFHs = [StarFormationHistories(s, fileorigin) for s in sfhs]

    # write trilegal sfh
    tri_sfhs = [S.make_trilegal_sfh(**make_trilegal_sfh_kw) for S in SFHs]
    return tri_sfhs


def vary_sfhs_of_one_galaxy(galaxy_name, cmd_input_file, mk_tri_sfh_kw=None,
                            match_sfh_file='default', vary_sfh_kw=None,
                            match_fileorigin='match-grid', make_many_kw=None,
                            galaxy_input_file='default', clean_first=False,
                            add_stage_lfs=None, ast=False, table_file='default'):
    '''
    Run a number of SFH variations on a galaxy.
    If passed default to the args, will attempt to find the file based on the
    galaxy name.

    ARGS:
    galaxy_name: target name, ex: ddo71 (case doesn't matter)
    cmd_input_file: filename, ex: 'cmd_input_CAF09_S_OCT13.dat'
    match_sfh_file: 'default', the sfh file from match.
    match_fileorigin: 'match-grid', which type of sfh file 'match' or 'match-grid'
    galaxy_input_file: 'default', base input file for trilegal, will be copied.

    mk_tri_sfh_kw: A dict to be passed to VarySFH.make_trilegal_sfh
        default: random_sfh = True, random_z = False

    make_many_kw: A dict to be passed to VarySFH.prepare_trilegal_sfr
        default: nsfhs = 50, mk_tri_sfh_kw dict.

    vary_sfh_kw: A dict to be passed to VarySFH.vary_the_sfh
        default: diag_plots = True, make_many_kw dict

    RETURNS:
    VarySFHs class
    '''
    # load and possibly overwrite make_many_kw from defaults
    make_many_kw = make_many_kw or {}
    vary_sfh_kw = vary_sfh_kw or {}
    mk_tri_sfh_kw = mk_tri_sfh_kw or {}

    mk_tri_sfh_kw = dict({'random_sfr': True, 'random_z': False}.items()
        + mk_tri_sfh_kw.items())

    make_many_kw = dict({'nsfhs': 50, 'mk_tri_sfh_kw': mk_tri_sfh_kw}.items()
        + make_many_kw.items())

    vary_sfh_kw = dict({'diag_plots': True, 'clean_first': clean_first,
                        'make_many_kw': make_many_kw}.items() + vary_sfh_kw.items())

    # load input files if not supplied
    if match_sfh_file == 'default':
        match_sfh_src = snap_src + '/data/sfh_parsec/match_files/'
        match_sfh_file, = rsp.fileIO.get_files(match_sfh_src, '%s*sfh' % galaxy_name.lower())

    if galaxy_input_file == 'default':
        galaxy_input_src = snap_src + '/input/'
        galaxy_input_file, = rsp.fileIO.get_files(galaxy_input_src, '*%s*dat' % galaxy_name.upper())

    # initialize VarySFHs object, loads ASTs, obs data
    vSFH = VarySFHs(galaxy_input_file, match_sfh_file, match_fileorigin,
                    target=galaxy_name.lower(), ast=ast, table_file=table_file)

    # vary the SFH
    vSFH.vary_the_SFH(cmd_input_file, add_stage_lfs=add_stage_lfs, **vary_sfh_kw)

    return vSFH


def simulation_from_beginning():
    import pdb
    pdb.set_trace()
    ags = AncientGalaxies()
    galaxy_table = ags.write_trgb_table(exclude_region=[0.1, 0.2])
    dry_run = False
    clean_first = True
    if dry_run is True:
        clean_first = False
    #targets = galaxy_tests.ancients()
    targets = ['ddo71']
    hist_it_up = False
    cmd_inputs = ['cmd_input_CAF09_S_NOV13.dat', 'cmd_input_CAF09_S_OCT13.dat']
    for cmd_input in cmd_inputs:
        for target in targets:
            vary_sfhs_of_one_galaxy(target, cmd_input, clean_first=clean_first,
                                    make_many_kw={'nsfhs': 50},
                                    vary_sfh_kw={'dry_run': dry_run,
                                                 'hist_it_up': hist_it_up},
                                    add_stage_lfs='default',
                                    table_file=galaxy_table)

if __name__ == '__main__':
    simulation_from_beginning()
