import os
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import ImageGrid
import random
import logging
logger = logging.getLogger()
import time
import pprint

from astroML.stats import binned_statistic_2d
import brewer2mpl

import ResolvedStellarPops as rsp
import ResolvedStellarPops.convertz as convertz
from TPAGBparams import table_src, snap_src
import galaxy_tests

angst_data = rsp.angst_tables.AngstTables()


def find_contamination_by_phases(output_files=None):
    if output_files is None:
        output_files = ['/home/phil/research/TP-AGBcalib/SNAP/models/varysfh/ddo71/caf09_s_nov13/mc/output_ddo71_caf09_s_nov13.dat',
                        #'/home/phil/research/TP-AGBcalib/SNAP/models/varysfh/ddo78/caf09_s_nov13/mc/output_ddo78_caf09_s_nov13.dat',
                        '/home/phil/research/TP-AGBcalib/SNAP/models/varysfh/hs117/caf09_s_nov13/mc/output_hs117_caf09_s_nov13.dat',
                        '/home/phil/research/TP-AGBcalib/SNAP/models/varysfh/kdg73/caf09_s_nov13/mc/output_kdg73_caf09_s_nov13.dat',
                        '/home/phil/research/TP-AGBcalib/SNAP/models/varysfh/kkh37/caf09_s_nov13/mc/output_kkh37_caf09_s_nov13.dat',
                        '/home/phil/research/TP-AGBcalib/SNAP/models/varysfh/ngc2976-deep/caf09_s_nov13/mc/output_ngc2976-deep_caf09_s_nov13.dat',
                        '/home/phil/research/TP-AGBcalib/SNAP/models/varysfh/ngc404/caf09_s_nov13/mc/output_ngc404_caf09_s_nov13.dat']
    comp_tab = read_completeness_table()
    for output_file in output_files:
        target = output_file.split('output_')[1].split('_')[0]
        sgal = rsp.Galaxies.simgalaxy(output_file, filter1=get_filter1(target), filter2='F814W')
        sgal.target = target
        opt_trgb = angst_data.get_item(target.upper(), 'mTRGB', extra_key=['%s,F814W' % get_filter1(target)])
        ir_trgb = angst_data.snap_tab3['target' == target.upper()]['mTRGB_F160W']
        comp90s = [comp_tab['target' == sgal.target.upper()]['opt_filter2'],
                   comp_tab['target' == sgal.target.upper()]['ir_filter2']]
        bright_limit = [opt_trgb-0.1, ir_trgb-0.25]
        contamination_by_phases(sgal, dim_limit=comp90s, bright_limit=bright_limit)
        #exclude_low = [opt_trgb+.1, ir_trgb+.25]
        #contamination_by_phases(sgal, dim_limit=exclude_low, bright_limit=bright_limit)
    return sgal


def default_output_location(target, extra_directory=None, mc=False):
    if extra_directory is None:
        outfile_loc = os.path.join(snap_src, 'models', 'varysfh', target)
        extra_directory = ''
    else:
        outfile_loc = os.path.join(snap_src, 'models', 'varysfh', target,
                                   extra_directory)
        extra_directory += '_'
    if mc is True:
        outfile_loc += '/mc/'
    else:
        outfile_loc += '/'
    return outfile_loc


def default_agb_filepath(cmd_input_file, extra_directory='default'):
    agb_model = cmd_input_file.replace('cmd_input_', '').lower()
    if extra_directory == 'default':
        agb_mod = agb_model.split('.')[0]
    else:
        agb_mod = None
    return agb_mod


def completeness_table(targets, comp_val, ast_kw=None, outfile='default'):

    targets = galaxy_tests.load_targets(targets)
    if outfile == 'default':
        outfile = os.path.join(snap_src + '/tables/', 'completeness_%.2f.dat' % comp_val)

    fmt = '%(target)s %(opt_filter1)s %(opt_filter2)s %(ir_filter1)s %(ir_filter2)s \n'
    with open(outfile, 'w') as out:
        out.write('# completeness fraction: %.2f \n' % comp_val)
        out.write('# ' + fmt.replace('%(','').replace(')s','').replace('\'',''))
        for target in targets:
            ast_dict = find_completeness(target, comp_val)
            ast_dict['target'] = target
            out.write(fmt % ast_dict)
    logger.info('wrote %s' % outfile)


def find_completeness(target, comp_val, ast_kw=None):
    ast_kw = ast_kw or {}
    ast_kw = dict({'combined_filters': True,
                   'interpolate': True}.items() + ast_kw.items())
    ast_dict = {}

    fake_files = galaxy_tests.get_fake_files(target)
    for band, fake_file in zip(['opt', 'ir'], fake_files):
        asts = rsp.Galaxies.artificial_star_tests(fake_file)
        asts.completeness(**ast_kw)
        (ast_dict['%s_filter1' % band], ast_dict['%s_filter2' % band]) = \
            asts.get_completeness_fraction(comp_val)
    return ast_dict


def read_completeness_table(table='default'):
    if table == 'default':
        table = snap_src + '/tables/completeness_0.90.dat'
    dtype = [('target', '|S16'), ('opt_filter1', '<f8'), ('opt_filter2', '<f8'),
             ('ir_filter1', '<f8'), ('ir_filter2', '<f8')]
    data = np.genfromtxt(table, dtype=dtype)
    return data


def number_of_stars(gal=None, exclude_region='default', mag_below_trgb=2.,
                    galaxy_information=None, factor=2., comp_frac=False,
                    indices=False):
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
        mag1 = galaxy_information.get('mag1')
        filter2 = galaxy_information.get('filter2')
        filter1 = galaxy_information.get('filter1')
    else:
        target = gal.target
        mag2 = gal.mag2
        mag1 = gal.mag1
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

    if comp_frac is False:
        offset = trgb + mag_below_trgb
    else:
        offset = mag_below_trgb
    # rgb will go from mag_below_trgb to trgb + exclude_region
    color_cut, = np.nonzero((mag1-mag2) > 0.3)
    nrgb = rsp.math_utils.between(mag2[color_cut], offset,
                                      trgb + exclude_region)

    # agb will go from trgb - exclude_region to very bright
    nagb = rsp.math_utils.between(mag2[color_cut], trgb - exclude_region, -99.)

    # add trgb_err to galaxy instance, could also add factor...
    if gal is not None:
        gal.trgb_err = trgb_err

    exclude_dict = {'trgb': trgb, 'trgb_err': trgb_err, 'factor': factor}
    if indices is False:
        nrgb = len(nrgb)
        nagb = len(nagb)
        
    return nrgb, nagb, exclude_dict
    
class PHATFields(object):
    def __init__(self):
        pass

    def write_trgb_table(self, brick=21, field=6, mag_below_trgb=[2, 1.5],
                         exclude_region='default', comp_table='default'):
        pass


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
                         exclude_region='default', comp_table='default'):
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
            logger.error('need to write how to read other targets')

        comp_frac = False
        if mag_below_trgb == 'comp_frac':
            logger.info('using completeness limit for rgb norm mag.')
            comp_frac = True
            comp_data = read_completeness_table(table=comp_table)
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
                if comp_frac is True:
                    mag_below_trgb = rsp.fileIO.item_from_row(comp_data,
                                                              'target', gal.target,
                                                              '%s_filter2' % band)
                else:
                    mag_below_trgb = mag_below_trgb[i]

                nrgb, nagb, ex_dict = number_of_stars(gal, comp_frac=comp_frac,
                                                      mag_below_trgb=mag_below_trgb,
                                                      exclude_region=exreg)
                mag_max = gal.mag2.max()
                mag_min = gal.mag2.min()
                nbins = np.sqrt(len(gal.mag2))
                gdict = {'%s_trgb' % band: gal.trgb,
                         '%s_trgb_err' % band: ex_dict['trgb_err'],
                         'n%s_rgb' % band: nrgb,
                         'n%s_agb' % band: nagb,
                         'n%s_bins' % band: nbins,
                         '%s_max' % band: mag_max,
                         '%s_min' % band: mag_min}
                try:
                    gal_dict[gal.target].update(gdict)
                except KeyError:
                    gal_dict[gal.target] = gdict
            factors.append(ex_dict['factor'])
        fmt = '%(target)s '
        fmt += '%(opt_trgb).2f %(opt_trgb_err).2f %(nopt_rgb)i %(nopt_agb)i '
        fmt += '%(ir_trgb).2f %(ir_trgb_err).2f %(nir_rgb)i %(nir_agb)i '
        fmt += '%(opt_max).2f %(ir_max).2f %(opt_min).2f %(ir_min).2f %(nopt_bins)i %(nir_bins)i \n'
        outfile = os.path.join(table_src, '%s_galaxies.dat' % tstring)
        if comp_frac is False:
            header = '# mags below trgb: optical=%.2f nir=%.2f \n' % (mag_below_trgb[0],
                                                                      mag_below_trgb[1])
        else:
            header = '# mags below trgb set by completeness fraction. \n'
        if len(np.unique(factors)) > 1:
            header += '# excluded trgb region: %.2f * opt trgb mag_err \n' % factors[0]
            header += '# excluded trgb region: %.2f * ir trgb mag_err \n' % factors[1]
        else:
            header += '# excluded trgb region: %.2f * trgb mag_err \n' % factors[0]
        header += '# '+ fmt.replace('%(', '').replace(').2f', '').replace(')i', '').replace(')s', '')

        with open(outfile, 'w') as out:
            out.write(header)
            for k, v in gal_dict.items():
                gal_dict[k]['target'] = k
                out.write(fmt % gal_dict[k])
        logger.info('wrote %s' % outfile)
        return outfile

    def read_trgb_table(self, table_name):
        with open(table_name, 'r') as tf:
            col_keys = tf.readlines()[3].replace('#','').strip().split()
        dtype = [(c, '<f8') for c in col_keys]
        dtype[0] = ('target', '|S16')
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
        if 'completeness' in mline:
            self.offsets = []
        else:
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
    if 'match-grid' == file_origin.lower() or 'match-hmc' == file_origin.lower():
        data = rsp.match_utils.read_binned_sfh(filename)
    elif 'match-old' == file_origin.lower():
        data = rsp.match_utils.read_match_old(filename)
    else:
        logger.error('please add a new data reader')
    
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
        self.file_origin = file_origin

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
                logger.warning('no uncertainties')
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

        if self.file_origin == 'match-old':
            random_sfr = False
            random_z = False
            zdisp = False
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
        if not 'm' in attr_str:
            ax.errorbar(self.data.lagei, val_arr, [errm_arr, errp_arr],
                        linestyle='steps-mid', lw=2, color='r')
        ax.set_ylabel(ylab, fontsize=20)
        ax.set_xlabel('$\log {\\rm Age (yr)}$', fontsize=20)
        ax.set_xlim(8, 10.5)
        if outfile is not None:
            plt.savefig(outfile, dpi=150)
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
            fig.savefig(outfig, dpi=150)


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
            table_file = os.path.join(table_src, 'ancients_0.1_0.2_galaxies.dat')
        else:
            logger.info('reading from table %s' % table_file)
        # Reads in the data as well as the mag offsets and factor in the
        # exclude region.
        self.ags = AncientGalaxies()
        self.ags.read_trgb_table(table_file)

        if len(self.ags.offsets) == 0:
            self.comp_data = read_completeness_table()

        # load stars in data RGB and AGB region as well as TRGBs.
        # this is just a call to self.ags.data
        self.load_data_for_normalization()

        # should we use ast corrections?
        self.ast = ast
        if self.ast is False:
            logger.info('not using artificial star corrections')

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
                logger.info('wrote %s' % new_out)
            self.galaxy_inputs.append(new_out)

    def prepare_outfiles(self, outfile_loc='default', extra_directory=None,
                         clean_first=False, extra_str=''):
        '''
        prepare outfiles for vary_the_sfh.
        opt_lf, ir_lf
        opt_lf_noagb, ir_lf_noagb
        opt_mass_met, ir_mass_met
        narratio

        NOTE: all attrs saved are file objects.
        '''
        if outfile_loc == 'default':
            self.outfile_loc = default_output_location(self.target,
                                                       extra_directory=extra_directory,
                                                       mc=self.mc)
        rsp.fileIO.ensure_dir(self.outfile_loc)
        logger.debug('outfile_loc set: %s' % self.outfile_loc)

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
            name_fmt = '%s_%s_%s%s.dat' % (extra_directory, self.target, fname,
                                           extra_str)
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
        if not hasattr(self, 'opt_mag'):
            self.load_trilegal_data()

        if hist_it_up is True:
            # hist it up!
            opt_hist, opt_bins = galaxy_tests.hist_it_up(self.opt_mag, threash=5)
            ir_hist, ir_bins = galaxy_tests.hist_it_up(self.ir_mag, threash=5)
        else:
            #nbins = (np.max(opt_mag) - np.min(opt_mag)) / 0.1
            opt_hist, opt_bins = np.histogram(self.opt_mag, bins=self.opt_bins)

            #nbins = (np.max(ir_mag) - np.min(ir_mag)) / 0.1
            ir_hist, ir_bins = np.histogram(self.ir_mag, bins=self.ir_bins)

        # scale the simulated LF to match the data LF
        opt_hist = np.array(opt_hist, dtype=float) * opt_norm
        ir_hist = np.array(ir_hist, dtype=float) * ir_norm

        np.savetxt(opt_file, (opt_hist, opt_bins[1:]), fmt='%g')
        np.savetxt(ir_file, (ir_hist, ir_bins[1:]), fmt='%g')

        logger.info('wrote %s, %s' % (opt_file.name, ir_file.name))
        return

    def count_stars_from_hist(self, opt_hist, opt_bins, ir_hist, ir_bins):
        ratio_data = {}
        for i, (hist, bins, band) in enumerate(zip([opt_hist, ir_hist],
                                                 [opt_bins, ir_bins],
                                                 ['opt', 'ir'])):
            trgb = self.__getattribute__('%s_trgb' % band)
            trgb_err = self.__getattribute__('%s_trgb_err' % band)
            norm = self.__getattribute__('%s_offset' % band)
            irgb = rsp.math_utils.between(bins, norm,
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
        logger.info('wrote %s' % out_file.name)
        logger.debug(pprint.pformat(out_dict))

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
            
        logger.info('wrote %s, %s' % (self.opt_mass_met_file.name,
                                      self.ir_mass_met_file.name))

    def load_data_for_normalization(self):
        '''load the numbers of data rgb and agb stars from self.ags'''
        column_names = self.ags.data.dtype.names
        if '404' in self.target:
            target = self.target.replace('-deep', '')
        else:
            target = self.target
        row, = np.nonzero(self.ags.data['target'] == target)
        try:
            [self.__setattr__('%s' % c, self.ags.data[row]['%s' % c][0])
             for c in column_names if c != 'target']
        except IndexError:
            '%s not found in data table' % target
        if len(self.ags.offsets) == 0:
            self.ir_offset = rsp.fileIO.item_from_row(self.comp_data,
                                                      'target', self.target.upper(),
                                                      'ir_filter2')
            self.opt_offset = rsp.fileIO.item_from_row(self.comp_data,
                                                      'target', self.target.upper(),
                                                      'opt_filter2')
        else:
            self.ir_offset = self.ir_trgb + self.ags.offsets[1]
            self.opt_offset = self.opt_trgb + self.ags.offsets[0]

        #self.opt_bins = np.linspace(self.opt_min, self.opt_max, self.nopt_bins)
        #if np.diff(self.opt_bins)[0] < 0.1:
        self.opt_bins = np.arange(self.opt_min, self.opt_max, 0.1)
        #self.ir_bins = np.linspace(self.ir_min, self.ir_max, self.nir_bins)
        #if np.diff(self.opt_bins)[0] < 0.1:
        self.ir_bins = np.arange(self.ir_min, self.ir_max, 0.1)
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

    def shift_mags(self, by_stage=True, opt_inds=None, ir_inds=None):
        '''shift mags to they agree with opt trgb'''
        opt_mag = self.sgal.data.get_col('F814W')
        ir_mag = self.sgal.data.get_col('F160W')
        if opt_inds is None:
            opt_inds = np.arange(len(opt_mag))
        if ir_inds is None:
            ir_inds = np.arange(len(ir_mag))


        if by_stage is True:
            # Threshold is set at 100 rgb stars in a bin.
            rgb_thresh = 50

            self.sgal.all_stages('RGB')
            rgb_inds = np.intersect1d(self.sgal.irgb, opt_inds)
            rgb_bins = np.arange(10, 30, 0.1)
            rgb_hist, _ = np.histogram(opt_mag[rgb_inds], bins=rgb_bins)
            rgb_bin_edge = np.nonzero(rgb_hist > rgb_thresh)[0][0] - 1
            opt_offset = rgb_bins[rgb_bin_edge] - self.opt_trgb

            rgb_inds = np.intersect1d(self.sgal.irgb, ir_inds)
            rgb_hist, _ = np.histogram(ir_mag[rgb_inds], bins=rgb_bins)
            rgb_bin_edge = np.nonzero(rgb_hist > rgb_thresh)[0][0] - 1
            ir_offset = rgb_bins[rgb_bin_edge] - self.ir_trgb

        else:
            # closest star in the trilegal catalog to star on the trgb
            trgb_color = angst_data.get_item(self.target, 'mean_color',
                                             extra_key=['%s,%s' % (self.filter1,
                                                                   'F814W')])
            ind, dist = rsp.math_utils.min_dist2d(trgb_color, self.opt_trgb,
                                                  opt_color, opt_mag)

            # correction for mag
            ir_offset = self.ir_trgb - ir_mag[ind]
            opt_offset = self.opt_trgb - opt_mag[ind]
        logger.debug('IR OFFSET: %f' % ir_offset)
        logger.debug('OPT OFFSET: %f' % opt_offset)

        self.opt_mag = opt_mag[opt_inds] - opt_offset
        self.ir_mag = ir_mag[ir_inds] - ir_offset

        return

    def mass_cut_inds(self):
        if not hasattr(self, 'mass_cut'):
            self.mass_cut = None
        else:
            logger.info('mass cut: %s' % self.mass_cut)
        mass = self.sgal.data.get_col('m_ini')
        if self.mass_cut is not None:
            nmasses = len(self.mass_cut.split('f')) - 1
            if nmasses == 1:
                inds, = np.nonzero([eval(self.mass_cut % m) for m in mass])
            if nmasses == 2:
                inds, = np.nonzero([eval(self.mass_cut % (m, m)) for m in mass])
        else:
            inds = np.arange(len(mass))
        return inds

    def load_trilegal_data(self):
        '''load trilegal F814W and F160W mags'''
        if not hasattr(self, 'ast'):
            self.ast = False

        extra = ''
        if self.ast is True:
            extra = '_cor'
        #print extra, 'extra'
        opt_mag = self.sgal.data.get_col('%s%s' % ('F814W', extra))
        ir_mag = self.sgal.data.get_col('%s%s' % ('F160W', extra))
        opt_mag1 = self.sgal.mag1
        ir_mag1 = self.sgal.data.get_col('%s%s' % ('F110W', extra))

        opt_color_cut, = np.nonzero((opt_mag1 - opt_mag) > 0.3)
        ir_color_cut, = np.nonzero((ir_mag1 - ir_mag) > 0.3)
        self.opt_color_cut = opt_color_cut
        self.ir_color_cut = ir_color_cut
        self.shift_mags(opt_inds=opt_color_cut, ir_inds=ir_color_cut)

        #inds = self.mass_cut_inds()

        # it's a matter of memory to not keep all the non-recovered stars
        # the hard coded limits just help not having to set ylims in LF plots.
        #if self.ast is True:
        #    opt_mag = opt_mag[np.isfinite(opt_mag)]
        #    opt_mag = opt_mag[opt_mag < 29]
        #    ir_mag = ir_mag[np.isfinite(ir_mag)]
        #    ir_mag = ir_mag[ir_mag < 26]
        return

def contamination_by_phases(self, sopt_rgb, sopt_agb, sir_rgb, sir_agb):
    regions = ['MS', 'RGB','RHEB', 'EAGB', 'TPAGB']
    self.sgal.all_stages()
    indss = [self.sgal.__getattribute__('i%s' % r.lower()) for r in regions]
    line = '%s ' % self.target
    print '#',' '.join(regions)
    
    for i, (rgb, agb, inds) in enumerate(zip([sopt_rgb, sir_rgb], [sopt_agb, sir_agb],
                                     [self.opt_color_cut, self.ir_color_cut])):
        #fig, (ax1, ax2) = plt.subplots(ncols=2)
        filter1 = self.filter1
        band = 'opt'
        if i == 1:
            filter1 = 'F160W'
            band = 'ir'
        ncontam_rgb = [np.intersect1d(inds[s], rgb) for s in indss]
        ncontam_agb = [np.intersect1d(inds[s], agb) for s in indss]
        
        #ax1.plot(color, sgal.data.get_col(filter2), '.')
        #ax2.plot(color, sgal.data.get_col(filter2), '.')
        #[ax1.plot(color[n], sgal.data.get_col(filter2)[n], '.',
        #         label=regions[j]) for j, n in enumerate(ncontam_rgb)]
        #[ax2.plot(color[n], mag2[n], '.') for n in ncontam_agb]       
        nrgb = np.max([1., len(ncontam_rgb[1])])
        line_rgb = 'rgb ' + band + ' ' + ' '.join(['%.3f' % (len(n))
                                                      for n in ncontam_rgb])
        nagb = np.max([1., len(ncontam_agb[-1])])
        
        line_agb = 'agb ' + band + ' ' + ' '.join(['%.3f' % (len(n))
                                                      for n in ncontam_agb])
        print line_rgb
        print line_agb
        #ax1.legend()
    print line

    def do_normalization(self, filter1=None, trilegal_output=None,
                         hist_it_up=True, dry_run=False):
        '''Do the normalization and save small part of outputs.'''
        if not hasattr(self, 'sgal'):
            assert trilegal_output is not None, \
                'need sgal loaded or pass trilegal catalog file name'
        if self.mc is True or not hasattr(self, 'sgal'):
            self.read_trilegal_catalog(trilegal_output, filter1=filter1)

        if not hasattr(self, 'opt_mag'):
            self.load_trilegal_data()

        # Recovered stars in simulated RGB region.

        sopt_rgb = self.sgal.stars_in_region(self.opt_mag,
                                             self.opt_offset,
                                             self.opt_trgb + \
                                             self.opt_trgb_err * \
                                             self.ags.factor[0])
        sir_rgb = self.sgal.stars_in_region(self.ir_mag,
                                            self.ir_offset,
                                            self.ir_trgb + \
                                            self.ir_trgb_err * self.ags.factor[1])

        # Recovered stars in simulated AGB region.
        sopt_agb = self.sgal.stars_in_region(self.opt_mag, self.opt_trgb - \
                                             self.opt_trgb_err * \
                                             self.ags.factor[0], 10.)
        sir_agb = self.sgal.stars_in_region(self.ir_mag, self.ir_trgb - \
                                            self.ir_trgb_err  * \
                                            self.ags.factor[1], 10.)

        # normalization
        opt_norm = self.nopt_rgb / float(len(sopt_rgb))
        ir_norm = self.nir_rgb / float(len(sir_rgb))
        logger.info('OPT Normalization: %f' % opt_norm)
        logger.info('IR Normalization: %f' % ir_norm)

        # for tracking non agb stars
        itpagb = self.sgal.stage_inds('TPAGB')
        non_tpagb = list(set(np.arange(len(self.sgal.mag2))) - set(itpagb))

        # random sample the data distribution
        rands = np.random.random(len(self.opt_mag))
        opt_ind, = np.nonzero(rands < opt_norm)
        rands = np.random.random(len(self.ir_mag))
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

        self.binary_contamination(opt_agb, ir_agb)

        #save LF in both filters
        #if self.mc is True:
        # only need to save the LF and the ratios.
        if dry_run is False:
            self.write_mass_met_file()
            self.write_LF(opt_norm, ir_norm, self.opt_lf_file, self.ir_lf_file,
                          hist_it_up=hist_it_up)
            self.write_LF(opt_norm, ir_norm, self.opt_lf_noagb_file, self.ir_lf_noagb_file,
                          hist_it_up=hist_it_up, inds=non_tpagb)
            self.write_ratio(opt_rgb, opt_agb, ir_rgb, ir_agb,
                             self.narratio_file)

        return opt_norm, ir_norm, opt_ind, ir_ind

    def binary_contamination(self, sopt_agb, sir_agb):
        '''fraction of binaries in trgb'''
        binaries, = np.nonzero(self.sgal.data.get_col('m2/m1') > 0)
        self.sgal.all_stages()
        rgb_binaries = np.intersect1d(self.sgal.irgb, binaries)
        nopt_binaries_above_trgb = float(len(np.intersect1d(sopt_agb, rgb_binaries)))
        nir_binaries_above_trgb = float(len(np.intersect1d(sir_agb, rgb_binaries)))

        tpagb_opt_contamination = nopt_binaries_above_trgb/float(len(sopt_agb))
        tpagb_ir_contamination = nir_binaries_above_trgb/float(len(sir_agb))

        logger.info('number of binaries: %i' % len(binaries))
        logger.info('number of rgb binaries in tpagb region opt, ir: %i, %i' % \
                    (nopt_binaries_above_trgb, nir_binaries_above_trgb))
        logger.info('frac rgb binaries in tpagb, region: %.3f, %.3f' % \
                    (tpagb_opt_contamination, tpagb_ir_contamination))

        return tpagb_opt_contamination, tpagb_ir_contamination

    def vary_the_SFH(self, cmd_input_file, make_many_kw=None, dry_run=False,
                     diag_plots=False, hist_it_up=True, outfile_loc='default',
                     extra_directory='default', clean_first=False,
                     add_stage_lfs=None, filter1=None, mass_cut=None,
                     extra_str=''):
        '''
        make the sfhs, make the galaxy inputs, run trilegal. For no trilegal
        runs, set dry_run True.
        '''
        self.mass_cut = mass_cut
        self.agb_mod = default_agb_filepath(cmd_input_file)

        self.mc = False
        nsfhs = make_many_kw.get('nsfhs', 50)
        if nsfhs > 1:
            self.mc = True

        self.prepare_outfiles(extra_directory=self.agb_mod,
                              clean_first=clean_first,
                              outfile_loc=outfile_loc,
                              extra_str=extra_str)

        make_many_kw = make_many_kw or {}
        if not 'mk_tri_sfh_kw' in make_many_kw.keys():
            make_many_kw['mk_tri_sfh_kw'] = {}

        make_many_kw['mk_tri_sfh_kw'] = dict({'outfile_fmt':
            self.sfr_outfilefmt, 'dry_run': dry_run}.items() +
            make_many_kw['mk_tri_sfh_kw'].items())

        # what's up to the logger
        add_file_logger(self.outfile_loc)
        logger.debug(pprint.pformat(locals))

        self.prepare_trilegal_sfr(make_many_kw=make_many_kw)

        self.prepare_galaxy_input(dry_run=dry_run)

        tname = 'output_%s_%s.dat' % (self.target, self.agb_mod)
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
                                  filter1=filter1, hist_it_up=hist_it_up,
                                  dry_run=dry_run)
        if dry_run is True:
            self.read_trilegal_catalog(trilegal_output, filter1=filter1)

        self.close_files()
        opt_chi2, ir_chi2 = self.poission_chi2(hist_it_up=hist_it_up)
        chi2_file = os.path.join(self.outfile_loc, '%s_chi2.dat' % self.target)
        with open(chi2_file, 'w') as c2:
            c2.write('# sfr opt_chi2 ir_chi2 \n')
            for i in range(len(opt_chi2)):
                c2.write('%s %.3f %.3f \n' % (self.sfr_files[i], opt_chi2[i], ir_chi2[i]))
        logger.info('wrote %s' % chi2_file)
        if diag_plots is True:
            outfile_fmt = os.path.join(self.outfile_loc, '%s_random%s' % (self.target, extra_str))
            outfile_fmt += '_%s.png'
            if make_many_kw['mk_tri_sfh_kw']['random_sfr'] is True:
                [self.plot_random_arrays(attr_str, from_files=True,
                                         outfile=outfile_fmt % attr_str)
                 for attr_str in ['sfr', 'mh']]

            self.compare_to_gal(hist_it_up=hist_it_up, add_stage_lfs=add_stage_lfs,
                                extra_str=extra_str)
            if self.mc is True:
                self.plot_mass_met_table(self.opt_mass_met_file.name,
                                         self.ir_mass_met_file.name)
    
    def plot_lf_file(self, opt_lf_file, ir_lf_file, axs=None, plt_kw=None,
                     opt_limit=None, ir_limit=None):
        '''needs work, but: plot the lf files.'''
        plt_kw = plt_kw or {}
        plt_kw = dict({'linestyle': 'steps',
                       'color': 'black',
                       'alpha': 0.2}.items() + plt_kw.items())
        if axs is None:
            fig, (axs) = plt.subplots(ncols=2, figsize=(12,6))
            plt.subplots_adjust(right=0.95, left=0.05, wspace=0.1)

        for i, (lf_file, limit) in enumerate(zip([opt_lf_file, ir_lf_file],
                                           [opt_limit, ir_limit])):

            hists, binss = self.load_lf_file(lf_file)

            for hist, bins in zip(hists, binss):
                if len(hist) != len(bins):
                    continue
                if limit is not None:
                    inds, = np.nonzero(bins <= limit)
                    axs[i].plot(bins[inds], hist[inds], **plt_kw)
                else:
                    axs[i].plot(bins, hist, **plt_kw)
        return axs

    def load_lf_file(self, lf_file):
        with open(lf_file, 'r') as lff:
            lines = [l.strip() for l in lff.readlines()
                     if not l.startswith('#')]

        hists = [np.array(l.split(), dtype=float) for l in lines[0::2]]
        binss = [np.array(l.split(), dtype=float) for l in lines[1::2]]
        return hists, binss

    def load_galaxies(self, hist_it_up=True):
        ir_gal = galaxy_tests.load_galaxy(self.target, band='ir')
        fits_src = snap_src + '/data/angst_no_trim'
        opt_gal = galaxy_tests.load_galaxy(self.target, band='opt',
                                           fits_src=fits_src)
        # make galaxy histograms
        if hist_it_up is True:
            opt_gal.hist, opt_gal.bins = galaxy_tests.hist_it_up(opt_gal.mag2,
                                                        threash=5)
            ir_gal.hist, ir_gal.bins = galaxy_tests.hist_it_up(ir_gal.mag2, threash=5)
        else:
            #nbins = (np.max(opt_gal.mag2) - np.min(opt_gal.mag2)) / 0.1
            opt_gal.hist, opt_gal.bins = np.histogram(opt_gal.mag2, bins=self.opt_bins)

            #nbins = (np.max(ir_gal.mag2) - np.min(ir_gal.mag2)) / 0.1
            ir_gal.hist, ir_gal.bins = np.histogram(ir_gal.mag2, bins=self.ir_bins)

        return opt_gal, ir_gal

    def compare_to_gal(self, hist_it_up=True, narratio=True, no_agb=False,
                       add_stage_lfs=None, extra_str=''):
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
        if narratio is True:
            # load ast_table
            ast_table = read_completeness_table()
        # load galaxy data
        opt_gal, ir_gal = self.load_galaxies(hist_it_up=hist_it_up)

        if self.mc is True:
            # plot lfs from simulations (and initialize figure)
            ax1, ax2 = self.plot_lf_file(self.opt_lf_file.name,
                                         self.ir_lf_file.name,
                                         opt_limit=opt_gal.comp50mag2,
                                         ir_limit=ir_gal.comp50mag2)
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
            if not hasattr(self, 'opt_mag'):
                self.load_trilegal_data()
            #stage_inds = self.mass_cut_inds()
            for ax, mag, norm, sinds in zip([ax1, ax2],
                                            [self.opt_mag, self.ir_mag],
                                            [self.opt_norm, self.ir_norm],
                                            [self.opt_color_cut, self.ir_color_cut]):
                self.sgal.make_lf(mag, stages=add_stage_lfs,
                                  hist_it_up=hist_it_up, stage_inds=sinds)
                for i in range(nstages):
                    try:
                        hist = self.sgal.__getattribute__('i%s_lfhist' % add_stage_lfs[i].lower())
                    except AttributeError:
                        continue
                    bins = self.sgal.__getattribute__('i%s_lfbins' % add_stage_lfs[i].lower())
                    stage_lf_kw['color'] = cols[i]
                    stage_lf_kw['label'] = '$%s$' % add_stage_lfs[i]
                    ax.plot(bins[:-1], hist*norm, **stage_lf_kw)

            sopt_hist, sopt_bins = self.sgal.make_lf(self.opt_mag, bins=self.opt_bins,
                                                     hist_it_up=hist_it_up)
            sir_hist, sir_bins = self.sgal.make_lf(self.ir_mag, bins=self.ir_bins,
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
        dplot_kw = {'linestyle': 'steps', 'color': 'darkred', 'lw': 2,
                    'label': '$%s$' % self.target.upper()}
        ax2.plot(ir_gal.bins[:-1], ir_gal.hist, **dplot_kw)
        ax1.plot(opt_gal.bins[:-1], opt_gal.hist, **dplot_kw)

        # initialize add numbers to the plot
        if narratio is True:
            if add_stage_lfs is None or self.mc is True:
                # count stars from the saved file (which counts random asts)
                ratio_data = rsp.fileIO.readfile(self.narratio_file.name)
            else:
                # count stars from the loaded histograms
                ratio_data = self.count_stars_from_hist(sopt_hist, sopt_bins,
                                                        sir_hist, sir_bins)
            stext_kw = {'color': 'black', 'fontsize': 14, 'ha': 'center'}
            dtext_kw = {'color': dplot_kw['color'], 'fontsize': 14,
                        'ha': 'center'}

        for i, (ax, gal, offset, trgb_err, band) in enumerate(zip([ax1, ax2],
                                                   [opt_gal, ir_gal],
                                                   [self.opt_offset,
                                                    self.ir_offset],
                                                   [self.opt_trgb_err,
                                                   self.ir_trgb_err],
                                                   ['opt', 'ir'])):
            ax.set_yscale('log')
            ax.set_ylim(1, ax.get_ylim()[1])
            # set the max to be the brightest of the data ... not great.
            xmax = np.round(self.sgal.data.get_col('%s' % gal.filter2).min())
            ax.set_xlim(xmax, gal.comp50mag2)

            yarr = np.linspace(*ax.get_ylim())
            # vertical lines around the trgb exclude region
            ax.fill_betweenx(yarr, gal.trgb - trgb_err * self.ags.factor[i],
                             gal.trgb + trgb_err * self.ags.factor[i],
                             color='black', alpha=0.2)

            ax.vlines(gal.trgb, *ax.get_ylim(), color='black',
                      linestyle='--')
            # line at dim limit for rgb normalization
            ax.vlines(offset, *ax.get_ylim(), linestyle='--',
                      color='black')
            
            # % completeness limit
            ast_frac = rsp.fileIO.item_from_row(ast_table, 'target',
                                                self.target.upper(),
                                                '%s_filter2' % band)
            ax.fill_betweenx(yarr, ast_frac, ax.get_xlim()[1],
                             color='black', alpha=0.2)
            ax.legend(loc=4, frameon=False)
            ax.set_xlabel('$%s$' % gal.filter2, fontsize=20)

            if narratio is True:
                # need to load the data nrgb and nagb, calculate the ratio
                # and error.
                nrgb = rsp.fileIO.item_from_row(self.ags.data, 'target',
                                                self.target, 'n%s_rgb' % band)
                nagb = rsp.fileIO.item_from_row(self.ags.data, 'target',
                                                self.target, 'n%s_agb' % band)
                dratio = nagb / nrgb
                dratio_err = galaxy_tests.count_uncert_ratio(nagb, nrgb)

                #yval = 1.2  # text yloc found by eye, depends on fontsize
                stext_kw['transform'] = ax.transAxes
                dtext_kw['transform'] = ax.transAxes
                yval = 0.95
                xagb_val = 0.17
                xrgb_val = 0.5
                xratio_val = 0.83
                xvals = [xagb_val, xrgb_val, xratio_val]

                # simulated nrgb and nagb are the mean values
                srgb_text = '$\langle N_{\\rm RGB}\\rangle =%i$' % \
                            np.mean(ratio_data['n%s_rgb' % band])
                sagb_text = '$\langle N_{\\rm TP-AGB}\\rangle=%i$' % \
                            np.mean(ratio_data['n%s_agb' % band])

                # one could argue taking the mean isn't the best idea for
                # the ratio errors.
                sratio_text = '$f=%.3f\pm%.3f$' % \
                              (np.mean(ratio_data['%s_ar_ratio' % band]),
                               np.mean(ratio_data['%s_ar_ratio_err' % band]))

                drgb_text = '$N_{\\rm RGB}=%i$' % nrgb
                dagb_text = '$N_{\\rm TP-AGB}=%i$' % nagb
                dratio_text =  '$f = %.3f\pm%.3f$' % (dratio, dratio_err)

                textss = [[srgb_text, sagb_text, sratio_text],
                         [drgb_text, dagb_text, dratio_text]]
                kws = [stext_kw, dtext_kw]

                for kw, texts in zip(kws, textss):
                    for xval, text in zip(xvals, texts):
                        ax.text(xval, yval, text, **kw)
                    yval -= .05  # stack the text

        ax1.set_ylabel('$\#$', fontsize=20)
        plt.tick_params(labelsize=16)
        outfile = os.path.join(self.outfile_loc, '%s_lfs%s.png' % (self.target, extra_str))
        plt.savefig(outfile, dpi=150)
        logger.info('wrote %s' % outfile)
        return ax1, ax2

    def plot_mass_met_table(self, opt_mass_met_file, ir_mass_met_file,
                            extra_str=''):
        fig = plt.figure(figsize=(8, 8))
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
                                '%s_mass_met%s.png' % (self.target, extra_str)), dpi=150)
        return grid

    def remove_files(self):
        logger.info('backing up everything in %s first' % self.outfile_loc)
        files = rsp.fileIO.get_files(self.outfile_loc, '*.*')
        bkdir = os.path.join(self.outfile_loc, 'bkup/')
        if os.path.isdir(bkdir):
            old_files = os.listdir(bkdir)
            if len(old_files) > 0:
                logger.info('overwriting a previous back up.')
        else:
            rsp.fileIO.ensure_dir(bkdir)
        [os.system('mv %s %s' % (f, bkdir)) for f in files]

    def poission_chi2(self, hist_it_up=False, lf_file_directory=None):
        if not hasattr(self, 'ir_bins'):
            self.load_data_for_normalization()
        if not hasattr(self, 'opt_lf_file'):
            self.mc = True
            self.prepare_outfiles(extra_directory=lf_file_directory)
            self.close_files()

        opt_gal, ir_gal = self.load_galaxies(hist_it_up=hist_it_up)

        # cut LF at 90% completeness
        obins, = np.nonzero(opt_gal.bins <= self.opt_offset)
        ibins, = np.nonzero(ir_gal.bins <= self.ir_offset)

        opt_model_hists, opt_models_binss = self.load_lf_file(self.opt_lf_file.name)
        ir_model_hists, ir_models_binss = self.load_lf_file(self.ir_lf_file.name)
        opt_chi2 = np.array([])
        ir_chi2 = np.array([])
        for i in range(len(ir_model_hists)):
            chi2, pct_dif, sig = rsp.Galaxies.stellar_prob(opt_gal.hist[obins[1:]],
                                                           opt_model_hists[i][obins[1:]])
            opt_chi2 = np.append(opt_chi2, chi2)
            chi2, pct_dif, sig = rsp.Galaxies.stellar_prob(ir_gal.hist[ibins[1:]],
                                                           ir_model_hists[i][ibins[1:]])
            ir_chi2 = np.append(ir_chi2, chi2)
        return opt_chi2, ir_chi2


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
                            add_stage_lfs=None, ast=False, table_file='default',
                            outfile_loc='default'):
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

    vary_sfh_kw = dict({'diag_plots': False, 'clean_first': clean_first,
                        'make_many_kw': make_many_kw}.items() + vary_sfh_kw.items())
    # load input files if not supplied
    if match_sfh_file == 'default':
        match_sfh_src = snap_src + '/data/sfh_parsec/'
        match_sfh_file, = rsp.fileIO.get_files(match_sfh_src, '%s*sfh' % galaxy_name.lower())

    elif match_sfh_file == 'grid':
        match_sfh_src = snap_src + '/data/sfh_parsec/match_grid/match_files/'
        match_sfh_file, = rsp.fileIO.get_files(match_sfh_src, '%s*sfh' % galaxy_name.lower())

    if galaxy_input_file == 'default':
        galaxy_input_src = snap_src + '/input/'
        galaxy_input_file, = rsp.fileIO.get_files(galaxy_input_src, '*%s*dat' % galaxy_name.upper())

    logger.debug('galaxy_input_file set: %s' % galaxy_input_file)
    logger.debug('match_sfh_file set: %s' % match_sfh_file)

    # initialize VarySFHs object, loads ASTs, obs data
    vSFH = VarySFHs(galaxy_input_file, match_sfh_file, match_fileorigin,
                    target=galaxy_name.lower(), ast=ast, table_file=table_file)

    # vary the SFH
    vSFH.vary_the_SFH(cmd_input_file, outfile_loc=outfile_loc,
                      add_stage_lfs=add_stage_lfs, **vary_sfh_kw)

    return vSFH


class StatisticalComparisons(object):
    def __init__(self, cmd_input_file):
        self.agb_mod = default_agb_filepath(cmd_input_file)

    def chi2tests(self, targets):
        targets = galaxy_tests.load_targets(targets)
        chi_dict = {}
        for target in targets:
            data = self.load_chi2files(target.lower())
            if len(data) == 0:
                continue
            for band in ['opt', 'ir']:
                field = '%s_chi2' % band
                chi_dict['%s_%s_mean' % (target, band)] = np.mean(data[field])
                chi_dict['%s_%s_std' % (target, band)] = np.std(data[field])
                print '%s %s $%.2f\pm%.2f$' % (target.upper(), band, np.mean(data[field]), np.std(data[field]))
            opt_total = [v for k,v in chi_dict.items() if 'opt' in k and 'mean' in k]
            ir_total = [v for k,v in chi_dict.items() if 'ir' in k and 'mean' in k]
            chi_dict['%s_opt_mean' % self.agb_mod] = np.mean(opt_total)
            chi_dict['%s_ir_mean' % self.agb_mod] = np.mean(ir_total)
            chi_dict['%s_opt_std' % self.agb_mod] = np.std(opt_total)
            chi_dict['%s_ir_std' % self.agb_mod] = np.std(ir_total)

            print '%s opt $%.2f\pm%.2f$' % (self.agb_mod, np.mean(opt_total), np.std(opt_total))
            print '%s ir $%.2f\pm%.2f$' % (self.agb_mod, np.mean(ir_total), np.std(ir_total))
        return chi_dict

    def load_chi2files(self, target, filefmt='%s_chi2.dat'):
        self.outfile_loc = default_output_location(target,
                                                   extra_directory=self.agb_mod,
                                                   mc=True)
        dtype = [('sfr_file', '|S130'), ('opt_chi2', '<f8'), ('ir_chi2', '<f8')]
        fname = os.path.join(self.outfile_loc, filefmt % target)
        if not os.path.isfile(fname):
            print '%s not found' % fname
            return np.array([])
        data = np.genfromtxt(fname, dtype=dtype)
        return data

    def narratio(self, targets):
        narr_dict = {}
        targets = galaxy_tests.load_targets(targets)
        for target in targets:
            data = self.load_narratio_file(target)
            if len(data) == 0:
                continue
            for band in ['opt', 'ir']:
                ratio = data['n%s_agb' % band]/data['n%s_rgb' % band]
                narr_dict['%s_%s_mean_ratio' % (target, band)] = np.mean(ratio)
                narr_dict['%s_%s_std' % (target, band)] = np.std(ratio)
                print '%s %s $%.2f\pm%.2f$' % (target.upper(), band,
                                               np.mean(ratio), np.std(ratio))
            opt_total = [v for k,v in narr_dict.items() if 'opt' in k and 'mean' in k]
            ir_total = [v for k,v in narr_dict.items() if 'ir' in k and 'mean' in k]
            narr_dict['%s_opt_mean' % self.agb_mod] = np.mean(opt_total)
            narr_dict['%s_ir_mean' % self.agb_mod] = np.mean(ir_total)
            print '%s opt $%.2f\pm%.2f$' % (self.agb_mod, np.mean(opt_total), np.std(opt_total))
            print '%s ir $%.2f\pm%.2f$' % (self.agb_mod, np.mean(ir_total), np.std(ir_total))

        return narr_dict

    def load_narratio_file(self, target):
        self.outfile_loc = default_output_location(target.lower(),
                                                   extra_directory=self.agb_mod,
                                                   mc=True)
        filefmt = '%s_%s_narratio.dat' % (self.agb_mod, target.lower())
        fname = os.path.join(self.outfile_loc, filefmt)
        if not os.path.isfile(fname):
            print '%s not found' % fname
            return np.array([])
        data = rsp.fileIO.readfile(fname)
        return data


def get_data(table_file='default'):
    if table_file == 'default':
        table_file = snap_src + '/tables/ancients_0.1_0.2_galaxies.dat'
    ags = AncientGalaxies()
    ags.read_trgb_table(table_file)
    data_dict = {}
    for i, target in enumerate(ags.data.target):
        for band in ['opt', 'ir']:
            ratio = ags.data[i]['n%s_agb' % band] / \
                    ags.data[i]['n%s_rgb' % band]
            data_dict['%s_%s' % (target, band)] = ratio
            data_dict['%s_%s_err' % (target, band)] = \
                galaxy_tests.count_uncert_ratio(ags.data[i]['n%s_agb' % band],
                                                ags.data[i]['n%s_rgb' % band])
    return data_dict


def narratio_table(targets, cmd_input_files, table_file='default'):

    scs = [StatisticalComparisons(c) for c in cmd_input_files]
    data_dict = get_data(table_file=table_file)
    model_dicts = [scs[i].narratio(targets) for i in range(len(cmd_input_files))]
    fmt = '$%.3f\\pm%.3f$ & '
    targets = galaxy_tests.load_targets(targets)
    for band in ['opt', 'ir']:
        print band
        header = 'Target & $\\frac{N_{\\rmTP-AGB}}{N_{\\rm RGB}}$ Data & '
        for i in range(len(cmd_input_files)):
            header += '$\\frac{N_{\\rm TP-AGB}}{N_{\\rm RGB}}$ %s & ' % \
                default_agb_filepath(cmd_input_files[i]).replace('caf09_s_', '')
            header += 'Frac. Difference & '
        header = header[:-2] + '\\\\ \n \\hline \n'
        print header
        for target in targets:
            line = '%s &' % target.upper()
            if '404' in target:
                continue
            dnarr = data_dict['%s_%s' % (target.lower(), band)]
            derr = data_dict['%s_%s_err' % (target.lower(), band)]
            dstr = fmt % (dnarr, derr)
            line += dstr
            for model_dict in model_dicts:
                mnarr = model_dict['%s_%s_mean_ratio' % (target.upper(), band)]
                mnerr = model_dict['%s_%s_std' % (target.upper(), band)]
                mstr = fmt % (mnarr, mnerr)
                line += mstr
                pct_diff = (mnarr - dnarr) / dnarr
                pct_diff_err = pct_diff * (mnerr/mnarr + derr/dnarr)
                pdstr = fmt % (pct_diff, pct_diff_err)
                line += pdstr
            line = line[:-2] + '\\\\ \n'
            print line


def chi2_table(targets, cmd_input_files, table_file='default'):
    scs = [StatisticalComparisons(c) for c in cmd_input_files]
    model_dicts = [scs[i].chi2tests(targets) for i in range(len(cmd_input_files))]
    fmt = '$%.3f\\pm%.3f$ & '
    targets = galaxy_tests.load_targets(targets)
    bands = ['opt', 'ir']
    header = 'Target & '
    for cmd_input_file in cmd_input_files:
        for band in bands:
            agb_mod = default_agb_filepath(cmd_input_file).replace('caf09_s_', '')
            header += '\\chi^2 %s %s & ' % (agb_mod, band)
    header = header[:-2] + '\\\\ \n'
    print header
    for target in targets:
        line = '%s &' % target.upper()
        if '404' in target:
            continue
        for model_dict in model_dicts:
            for band in bands:
                line += fmt % (model_dict['%s_%s_mean' % (target.upper(), band)],
                               model_dict['%s_%s_std' % (target.upper(), band)])
        print line
    line = 'Mean $\\chi^2$ & '
    for i, model_dict in enumerate(model_dicts):
        agb_mod = default_agb_filepath(cmd_input_files[i])
        for band in bands:
            line += fmt % (model_dict['%s_%s_mean' % (agb_mod, band)],
                           model_dict['%s_%s_std' % (agb_mod, band)])
    print line


def add_file_logger(logdir):
    # setup logger
    logger = logging.getLogger()

    #logging.basicConfig(filename=logfile, level=logging.DEBUG)
    # file handler
    rsp.fileIO.ensure_dir(logdir)
    logfname = time.strftime("log_%Y%m%d_%H%M%S.dat", time.localtime())
    logfile = os.path.join(os.path.abspath(logdir), logfname)
    fh = logging.FileHandler(logfile)
    fh.setLevel(logging.DEBUG)

    # create formatter and add it to the handlers
    formatter = logging.Formatter('%(asctime)-15s %(levelname)s %(funcName)s %(lineno)d %(message)s')
    fh.setFormatter(formatter)

    # add the handlers to the logger
    logger.addHandler(fh)
    logger.info('logfile set: %s' % logfile)


def data_table(targets):
    # target, av, dist, [opt: frac comp, trgb, nrgb, nagb ratio] [ir: ...]
    table_file = snap_src + '/tables/ancients_0.1_0.2_galaxies.dat'
    ags = AncientGalaxies()
    ags.read_trgb_table(table_file)
    data_dict = get_data(table_file=table_file)
    comp_data = read_completeness_table()

    for target in targets:
        if target == 'NGC2976-DEEP':
            extra_key='F606W,F814W'
        else:
            extra_key=None
        (Av, dmod) = [angst_data.get_item(target, i, extra_key=extra_key) for i in ['Av', 'dmod']]
        comp_row = rsp.fileIO.get_row(comp_data, 'target', target)
        nstars_row = rsp.fileIO.get_row(ags.data, 'target', target)
        sub_dict  = {}
        for (k,v) in data_dict.items():
            if '404' in target:
                target = target.lower().replace('-deep', '')
            if target in k or target.lower() in k:
                sub_dict[k] = v

        row = '%s & %.2f & %.2f & ' % (target, Av, dmod)
        row += '%(opt_filter2).2f & ' % comp_row
        row += '%(opt_trgb).2f & %(nopt_agb)i & %(nopt_rgb)i & ' % nstars_row
        row += '$%.3f\\pm%.3f$ & ' % (sub_dict['%s_opt' % target.lower()], sub_dict['%s_opt_err' % target.lower()])

        row += '%(ir_filter2).2f & ' % comp_row
        row += '%(ir_trgb).2f & %(nir_agb)i & %(nir_rgb)i & ' % nstars_row
        row += '$%.3f\\pm%.3f$ \\\\' % (sub_dict['%s_ir' % target.lower()], sub_dict['%s_ir_err' % target.lower()])
        print row


def paolas_tests():
    '''tests paola requested'''
    targets = ['ddo78', 'ddo71']
    cmd_inputs = ['cmd_input_CAF09_S_NOV13.dat']
    nsfhs = 1
    mk_tri_sfh_kw = {'random_sfr': False}
    match_sfh_file = 'grid'

    common_dict = {'match_sfh_file': match_sfh_file,
                   'mk_tri_sfh_kw': mk_tri_sfh_kw}
    ##1 No Dust
    #print '\n No Dust \n '
    #galaxy_input_filefmt = snap_src + '/input/tests/input_%s_nodust.dat'
    #extra_str = '_nodust'
    #simulation_from_beginning(targets, cmd_inputs, nsfhs, extra_str=extra_str,
    #                          galaxy_input_filefmt=galaxy_input_filefmt,
    #                          **common_dict)
    ##2 C Dust
    #print '\n C Dust \n '
    #galaxy_input_filefmt = snap_src + '/input/tests/input_%s_cdust.dat'
    #extra_str = '_cdust'
    #simulation_from_beginning(targets, cmd_inputs, nsfhs, extra_str=extra_str,
    #                          galaxy_input_filefmt=galaxy_input_filefmt,
    #                          **common_dict)
    #
    ##2 M Dust
    #print '\n M Dust \n '
    #galaxy_input_filefmt = snap_src + '/input/tests/input_%s_mdust.dat'
    #extra_str = '_mdust'
    #simulation_from_beginning(targets, cmd_inputs, nsfhs, extra_str=extra_str,
    #                          galaxy_input_filefmt=galaxy_input_filefmt,
    #                          **common_dict)
    #2 Aringer Dust
    print '\n Aringer Dust \n '
    galaxy_input_filefmt = snap_src + '/input/tests/input_%s_aringer.dat'
    extra_str = '_aringer'
    simulation_from_beginning(targets, cmd_inputs, nsfhs, extra_str=extra_str,
                              galaxy_input_filefmt=galaxy_input_filefmt,
                              **common_dict)
    #3 Big Exclude
    #print '\n Big Exclude \n '
    #ags = AncientGalaxies()
    #galaxy_table = ags.write_trgb_table(exclude_region=[0.2, 0.4],
    #                                    mag_below_trgb='comp_frac')
    #extra_str = '_bigexclude'
    #simulation_from_beginning(targets, cmd_inputs, nsfhs, extra_str=extra_str,
    #                          galaxy_table=galaxy_table,
    #                          **common_dict)

    # 4 Mass Cut
    #print '\n Mass Cut \n '
    #mass_cut = '(%f<1.) | (%f>1.2)'
    #extra_str = '_masscut'
    #simulation_from_beginning(targets, cmd_inputs, nsfhs, extra_str=extra_str,
    #                          mass_cut=mass_cut,
    #                          **common_dict)

    # Normal
    #print '\n Default \n '
    #simulation_from_beginning(targets, cmd_inputs, nsfhs,
    #                          **common_dict)

def simulation_from_beginning(targets, cmd_inputs, nsfhs, hist_it_up=False,
                              dry_run=False, galaxy_input_filefmt=None,
                              galaxy_table='default', mass_cut=None,
                              extra_str='', clean_first=False,
                              mk_tri_sfh_kw=None, match_sfh_file='default',
                              debug=False):

    # start the stream logger. The file logger is done target by target.
    global logger
    logger = logging.getLogger()
    logger.info('start of run: %s' % time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))
    # create formatter and add it to the handlers
    formatter = '%(asctime)-15s %(levelname)s %(funcName)s %(lineno)d %(message)s'
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.info('start of run: %s' % time.strftime("%a, %d %b %Y %H:%M:%S", time.localtime()))

    # cant do a dry run with no files.
    clean_first = True
    if dry_run is True:
        clean_first = False

    if debug is True:
        import pdb; pdb.set_trace()

    for target in targets:
        if galaxy_input_filefmt is None:
            galaxy_input_file = 'default'
        else:
            galaxy_input_file = galaxy_input_filefmt % target.upper()
        if '404' in target:
            target = 'ngc404'
        for cmd_input in cmd_inputs:
            vary_sfhs_of_one_galaxy(target, cmd_input, clean_first=clean_first,
                                    mk_tri_sfh_kw=mk_tri_sfh_kw,
                                    make_many_kw={'nsfhs': nsfhs},
                                    vary_sfh_kw={'dry_run': dry_run,
                                                 'hist_it_up': hist_it_up,
                                                 'mass_cut': mass_cut,
                                                 'extra_str': extra_str},
                                    add_stage_lfs='default',
                                    table_file=galaxy_table,
                                    galaxy_input_file=galaxy_input_file,
                                    match_sfh_file=match_sfh_file)

    #data_table(targets)
    #narratio_table(targets, cmd_inputs, table_file=galaxy_table)
    #chi2_table(targets)

if __name__ == '__main__':
    #paolas_tests()
    ags = AncientGalaxies()
    galaxy_table = ags.write_trgb_table(exclude_region=[0.1, 0.2],
                                        mag_below_trgb='comp_frac')
    #targets = ['ddo78', 'ddo71', 'hs117', 'kkh37', 'ngc2976-deep', 'ngc404-deep']
    targets = ['ngc2976-deep', 'ngc404-deep']
    
    mk_tri_sfh_kw = {'random_sfr': True, 'random_z': False}
    simulation_from_beginning(targets, ['cmd_input_CAF09_S_NOV13.dat',
                                        'cmd_input_CAF09_S_NOV13eta0.dat',
                                        'cmd_input_CAF09_S_OCT13.dat'],
                              50, mk_tri_sfh_kw=mk_tri_sfh_kw, debug=True)
