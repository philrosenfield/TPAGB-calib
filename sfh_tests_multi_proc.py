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
import itertools
import tables

from model_plots import load_plot_limits

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


def get_color_cut(filter1):
    '''
    see snap_src + tables/color_cuts.dat
    '''
    if filter1 == 'F606W':
        color_cut = 0.2
    elif filter1 == 'F475W':
        color_cut = 0.3
    elif filter1 == 'F110W':
        color_cut = 0.1
    else:
        color_cut = np.nan
        print 'warning, filter not found, no color cut'
    return color_cut


def get_filter1(target, fits_src=None):
    '''get optical filter1'''
    fits_src = fits_src or snap_src + '/data/angst_no_trim'
    fitsfile, = rsp.fileIO.get_files(fits_src, '*%s*.fits' % target.lower())
    filter1 = fitsfile.split(target)[1].split('_')[1].split('-')[0].upper()
    return filter1


def default_output_location(target, extra_directory=None, mc=False):
    outfile_loc = os.path.join(snap_src, 'models', 'varysfh', target)

    if extra_directory is not None:
        outfile_loc = os.path.join(outfile_loc, extra_directory)

    if mc is True:
        outfile_loc += '/mc'

    outfile_loc += '/'
    rsp.fileIO.ensure_dir(outfile_loc)

    return outfile_loc


def default_agb_filepath(cmd_input_file, extra_directory='default'):
    agb_model = cmd_input_file.replace('cmd_input_', '').lower()
    if extra_directory == 'default':
        agb_mod = agb_model.split('.')[0]
    else:
        agb_mod = None
    return agb_mod


def setup_files(cmd_input_file, target, outfile_loc='default', extra_str='',
                mc=True):

    agb_mod = default_agb_filepath(cmd_input_file)

    if outfile_loc == 'default':
        outfile_loc = \
            default_output_location(target, extra_directory=agb_mod, mc=mc)

    names = ['opt_lf', 'ir_lf', 'narratio', 'opt_mass_met', 'ir_mass_met',
             'contam']
    name_fmt = '%s_%s_%s%s.dat'

    fnames = [os.path.join(outfile_loc,
                           name_fmt % (agb_mod, target, f, extra_str))
              for f in names]

    return outfile_loc, fnames, agb_mod


def load_default_ancient_galaxies(table_file='default'):
    if table_file is 'default':
        table_file = os.path.join(table_src, 'ancients_0.1_0.2_galaxies.dat')
    else:
        logger.info('reading from table %s' % table_file)
    # Reads in the data as well as the mag offsets and factor in the
    # exclude region.
    ags = AncientGalaxies()
    ags.read_trgb_table(table_file)
    return ags


def combine_list_of_dictionaries(dlist):
    result_dict = {}
    for dic in dlist:
        for key in dic.keys():
            if not key in result_dict.keys():
                result_dict[key] = []
            result_dict[key].append(dic[key])

    return result_dict


def default_values_for_vsfh(target, cmd_input, vsfh_kw=None,
                            match_sfh_src=None, match_sfh_file='default',
                            galaxy_input_src='default',
                            galaxy_input_search_fmt='*%s*dat'):

    vsfh_kw = vsfh_kw or {}
    if match_sfh_file == 'default':
        match_sfh_src = snap_src + '/data/sfh_parsec/'
    elif match_sfh_file == 'grid':
        # oldest version!!
        match_sfh_src = snap_src + '/data/sfh_parsec/match_grid/match_files/'

    if galaxy_input_src == 'default':
        galaxy_input_src = snap_src + '/input/'
    galaxy_input, = rsp.fileIO.get_files(galaxy_input_src,
                                         galaxy_input_search_fmt % target.upper())

    vsfh_kw['sfh_file'], = rsp.fileIO.get_files(match_sfh_src,
                                                '%s*sfh' % target.replace('-deep', ''))
    vsfh_kw['galaxy_input'] = galaxy_input
    vsfh_kw['target'] = target
    vsfh_kw['cmd_input_file'] = cmd_input

    return vsfh_kw


def prepare_vsfh_run(targets, cmd_input_files, nsfhs, mk_tri_sfh_kw=None,
                     vary_sfh_kw=None, make_many_kw=None, vsfh_kw=None,
                     table_file='default', dry_run=False, default_kw=None):
    '''
    Run a number of SFH variations on a galaxy.
    If passed default to the args, will attempt to find the file based on the
    galaxy name.

    ARGS:
    galaxy_name: target name, ex: ddo71 (case doesn't matter)
    cmd_input_file: filename, ex: 'cmd_input_CAF09_S_OCT13.dat'
    match_sfh_file: 'default', the sfh file from match.
    match_fileorigin: 'match-grid', which type of sfh file
                        'match' or 'match-grid'
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
    # start the stream logger. The file logger is done target by target.
    global logger
    logger = logging.getLogger()
    logger.info('start of run: %s' % time.strftime("%a, %d %b %Y %H:%M:%S",
                                                   time.localtime()))
    # create formatter and add it to the handlers
    formatter = \
        '%(asctime)-15s %(levelname)s %(funcName)s %(lineno)d %(message)s'
    ch = logging.StreamHandler()
    ch.setLevel(logging.ERROR)
    ch.setFormatter(formatter)
    logger.addHandler(ch)
    logger.info('start of run: %s' % time.strftime("%a, %d %b %Y %H:%M:%S",
                                                   time.localtime()))

    make_many_kw = make_many_kw or {}
    vary_sfh_kw = vary_sfh_kw or {}
    mk_tri_sfh_kw = mk_tri_sfh_kw or {}
    default_kw = default_kw or {}

    mk_tri_sfh_kw = dict({'dry_run': dry_run}.items() + mk_tri_sfh_kw.items())

    make_many_kw = dict({'mk_tri_sfh_kw': mk_tri_sfh_kw}.items()
                         + make_many_kw.items())

    vary_sfh_kw = dict({'make_many_kw': make_many_kw}.items()
                        + vary_sfh_kw.items())

    vsfh_kw = vsfh_kw or {}
    vsfh_kw = dict({'file_origin': 'match-grid', 'table_file': table_file,
                    'outfile_loc': 'default', 'nsfhs': nsfhs}.items() \
                   + vsfh_kw.items())

    vsfh_kws = []
    vSFHs = []
    for target, cmd_input in itertools.product(targets, cmd_input_files):
        target = target.lower()
        #print target, cmd_input
        vsfh_kws.append(default_values_for_vsfh(target, cmd_input,
                                                vsfh_kw=vsfh_kw, **default_kw))
        vSFHs.append(VarySFHs(**vsfh_kw))
    return vSFHs, vary_sfh_kw


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
    color_cut, = np.nonzero((mag1-mag2) > get_color_cut(filter1))
    nrgb = rsp.math_utils.between(mag2[color_cut], offset,
                                      trgb + exclude_region)

    # agb will go from trgb - exclude_region to very bright
    nagb = rsp.math_utils.between(mag2[color_cut], trgb - exclude_region, -99.)

    # add trgb_err to galaxy instance, could also add factor..
    if gal is not None:
        gal.trgb_err = trgb_err

    exclude_dict = {'trgb': trgb, 'trgb_err': trgb_err, 'factor': factor}
    if indices is False:
        nrgb = len(nrgb)
        nagb = len(nagb)

    return nrgb, nagb, exclude_dict


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
            comp_data = tables.read_completeness_table(table=comp_table)
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


class StarFormationHistories(object):
    '''Make TRILEGAL star formation history files from MATCH'''
    def __init__(self, sfh_file, file_origin, sfr_files=None, sfr_file_loc=None,
                 sfr_file_search_fmt=None):
        self.base, self.name = os.path.split(sfh_file)
        self.data = tables.parse_sfh_data(sfh_file, file_origin)
        self.file_origin = file_origin
        self.sfr_files = sfr_files
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

    def plot_sfh(self, attr_str, ax=None, outfile=None, yscale='linear',
                 plot_random_arrays_kw=None, errorbar_kw=None):
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

        if ax is None:
            fig, axs = plt.subplots(figsize=(8, 8), nrows=2, sharex=True)

        for ax in axs:
            if not 'm' in attr_str:
                ax.errorbar(self.data.lagei, val_arr, [errm_arr, errp_arr],
                            **errorbar_kw)

            if len(plot_random_arrays_kw) > 0:
                # if loading the random arrays from files, need to give the
                # attribute to load.
                if plot_random_arrays_kw['from_files'] is True:
                    plot_random_arrays_kw['attr_str'] = attr_str
                self.plot_random_arrays(ax=ax, **plot_random_arrays_kw)
            ax.set_ylabel(ylab, fontsize=20)
            ax.set_xlim(8, 10.13)

        axs[1].set_xlabel('$\log {\\rm Age (yr)}$', fontsize=20)
        axs[0].set_yscale('log')
        axs[0].set_ylabel(ylab.replace('rm', 'rm\ \log\ '))
        fig.subplots_adjust(hspace=0.07)
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

    def make_many_trilegal_sfhs(self, nsfhs=100, mk_tri_sfh_kw=None):
        '''
        make nsfhs number of trilegal sfh input files.
        '''
        mk_tri_sfh_kw = mk_tri_sfh_kw or {}

        if nsfhs > 1:
            self.mc = True

        outfile_fmt = mk_tri_sfh_kw.get('outfile_fmt', 'default')

        # need to iterate outside of dict.. see below.
        del mk_tri_sfh_kw['outfile_fmt']

        # update any passed kw to make_trilegal_sfh with defaults.
        mk_tri_sfh_kw = dict({'random_sfr': True, 'random_z': False}.items()
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


class FileIO(object):
    def __init__(self):
        pass

    def check_target(self, target):
        if target is None:
            assert hasattr(self, 'target'), \
                'need to pass target or have attribute'
        else:
            self.target = target
        return

    def load_data_for_normalization(self, target=None, ags=None):
        self.check_target(target)
        '''load the numbers of data rgb and agb stars from self.ags'''
        if ags is None:
            ags = self.ags

        if len(ags.offsets) == 0:
            self.comp_data = tables.read_completeness_table()

        column_names = ags.data.dtype.names
        if '404' in self.target:
            target = self.target.replace('-deep', '')
        else:
            target = self.target
        row, = np.nonzero(ags.data['target'] == target)
        try:
            [self.__setattr__('%s' % c, ags.data[row]['%s' % c][0])
             for c in column_names if c != 'target']
        except IndexError:
            '%s not found in data table' % target
        if len(ags.offsets) == 0:
            self.ir_offset = rsp.fileIO.item_from_row(self.comp_data,
                                                      'target', target.upper(),
                                                      'ir_filter2')
            self.opt_offset = rsp.fileIO.item_from_row(self.comp_data,
                                                      'target', target.upper(),
                                                      'opt_filter2')
        else:
            self.ir_offset = self.ir_trgb + ags.offsets[1]
            self.opt_offset = self.opt_trgb + ags.offsets[0]

        self.opt_bins = np.arange(self.opt_min, self.opt_max, 0.1)
        self.ir_bins = np.arange(self.ir_min, self.ir_max, 0.1)

        return self.nopt_rgb, self.nir_rgb

    def load_lf_file(self, lf_file):
        with open(lf_file, 'r') as lff:
            lines = [l.strip() for l in lff.readlines()
                     if not l.startswith('#')]

        hists = [np.array(l.split(), dtype=float) for l in lines[0::2]]
        binss = [np.array(l.split(), dtype=float) for l in lines[1::2]]
        return hists, binss

    def load_galaxies(self, hist_it_up=True, target=None, ags=None,
                      color_cut=False):
        self.check_target(target)
        if not hasattr(self, 'opt_bins'):
            self.load_data_for_normalization(ags=ags)

        ir_gal = galaxy_tests.load_galaxy(self.target, band='ir')
        fits_src = snap_src + '/data/angst_no_trim'
        opt_gal = galaxy_tests.load_galaxy(self.target, band='opt',
                                           fits_src=fits_src)
        # make galaxy histograms
        if color_cut is True:
            filter1 = get_filter1(target)
            opt_color_inds = np.nonzero(opt_gal.color > get_color_cut(filter1))
            ir_color_inds = np.nonzero(ir_gal.color > get_color_cut('F110W'))
            opt_gal.color_cut = opt_color_inds
            ir_gal.color_cut = ir_color_inds
        if hist_it_up is True:
            opt_gal.hist, opt_gal.bins = galaxy_tests.hist_it_up(opt_gal.mag2,
                                                        threash=5)
            ir_gal.hist, ir_gal.bins = galaxy_tests.hist_it_up(ir_gal.mag2,
                                                               threash=5)
        else:
            #nbins = (np.max(opt_gal.mag2) - np.min(opt_gal.mag2)) / 0.1
            if color_cut is False:
                opt_gal.hist, opt_gal.bins = np.histogram(opt_gal.mag2,
                                                          bins=self.opt_bins)
                ir_gal.hist, ir_gal.bins = np.histogram(ir_gal.mag2,
                                                        bins=self.ir_bins)
            else:
                opt_gal.hist, opt_gal.bins = \
                    np.histogram(opt_gal.mag2[opt_color_inds],
                                 bins=self.opt_bins)
                ir_gal.hist, ir_gal.bins = \
                    np.histogram(ir_gal.mag2[ir_color_inds], bins=self.ir_bins)

        return opt_gal, ir_gal

    def load_trilegal_data(self):
        '''load trilegal F814W and F160W mags'''

        if hasattr(self, 'target'):
            filter1 = get_filter1(self.target)
        else:
            print 'help, I need filter1!!'
        opt_mag = self.sgal.data.get_col('F814W')
        ir_mag = self.sgal.data.get_col('F160W')
        opt_mag1 = self.sgal.mag1
        ir_mag1 = self.sgal.data.get_col('F110W')

        opt_color_cut, = \
            np.nonzero((opt_mag1 - opt_mag) > get_color_cut(filter1))
        ir_color_cut, = \
            np.nonzero((ir_mag1 - ir_mag) > get_color_cut('F110W'))
        self.opt_color_cut = opt_color_cut
        self.ir_color_cut = ir_color_cut
        self.shift_mags(opt_inds=opt_color_cut, ir_inds=ir_color_cut)
        return

    def read_trilegal_catalog(self, trilegal_output, filter1='F606W'):
        '''read the trilegal cat mag1 and mag2 are optical.'''
        self.sgal = rsp.Galaxies.simgalaxy(trilegal_output, filter1=filter1,
                                           filter2='F814W')
        return

    def shift_mags(self, opt_inds=None, ir_inds=None):
        '''shift mags to they agree with opt trgb'''
        opt_mag = self.sgal.data.get_col('F814W')
        ir_mag = self.sgal.data.get_col('F160W')
        if opt_inds is None:
            opt_inds = np.arange(len(opt_mag))

        if ir_inds is None:
            ir_inds = np.arange(len(ir_mag))

        # Threshold is set at 100 rgb stars in a bin.
        rgb_thresh = 5.

        self.sgal.all_stages('RGB')
        rgb_inds = np.intersect1d(self.sgal.irgb, opt_inds)
        rgb_bins = np.arange(10, 30, 0.01)
        rgb_hist, _ = np.histogram(opt_mag[rgb_inds], bins=rgb_bins)
        rgb_bin_edge = np.nonzero(rgb_hist > rgb_thresh)[0][0] - 1
        opt_offset = rgb_bins[rgb_bin_edge] - self.opt_trgb

        rgb_inds = np.intersect1d(self.sgal.irgb, ir_inds)
        rgb_hist, _ = np.histogram(ir_mag[rgb_inds], bins=rgb_bins)
        rgb_bin_edge = np.nonzero(rgb_hist > rgb_thresh)[0][0] - 1
        ir_offset = rgb_bins[rgb_bin_edge] - self.ir_trgb

        # HERE's A HACK FOR NO OFFSETS!!!
        #ir_offset = 0.
        #opt_offset = 0.

        logger.debug('IR OFFSET: %f' % ir_offset)
        logger.debug('OPT OFFSET: %f' % opt_offset)
        self.ir_moffset = ir_offset
        self.opt_moffset = opt_offset
        print ir_offset
        print opt_offset
        self.opt_mag = opt_mag[opt_inds] - opt_offset
        self.ir_mag = ir_mag[ir_inds] - ir_offset
        return


class VarySFHs(StarFormationHistories, AncientGalaxies, FileIO):
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
    def __init__(self, galaxy_input, sfh_file, file_origin, cmd_input_file,
                 nsfhs, outfile_loc='default', filter1=None, extra_str='',
                 target=None, table_file='default', just_once=False):
        '''
        galaxy_input is a template.
        '''
        # load SFH instance to make lots of trilegal runs
        if just_once is False:
            StarFormationHistories.__init__(self, sfh_file, file_origin)

        # add information to self, in other words, save
        # to use in plotting and stats later
        if target is None:
            gname = os.path.split(galaxy_input)[1]
            target = gname.split('_')[1].replace('.dat', '').lower()
        self.target = target
        self.filter1 = get_filter1(self.target)
        self.galaxy_input = galaxy_input
        self.cmd_input_file = cmd_input_file
        self.nsfhs = nsfhs
        self.mc = False
        if self.nsfhs > 1:
            self.mc = True
        self.extra_str = extra_str

        # exclude regions and the number of data rgb and agb stars
        self.ags = load_default_ancient_galaxies(table_file=table_file)
        self.load_data_for_normalization()

        # 90% (or whatever) completeness magnitudes
        if len(self.ags.offsets) == 0:
            self.comp_data = tables.read_completeness_table()

        # setup the locations all the files to write and read from
        self.outfile_loc, self.fnames, self.agb_mod = \
            setup_files(cmd_input_file, outfile_loc=outfile_loc,
                        extra_str=extra_str, mc=self.mc, target=self.target)

        # header files are needed in two cases
        # nagb/nrgb ratio file
        self.narratio_header = '# target nopt_rgb nopt_agb nir_rgb nir_agb '
        self.narratio_header += 'opt_ar_ratio ir_ar_ratio opt_ar_ratio_err '
        self.narratio_header += 'ir_ar_ratio_err \n'

        # contamination of phases in rgb and agb region file
        # (changing the self.regions will disrupt calculation of
        # rheb_eagb_contamination -- see contamination_by_phases code)
        self.regions = ['MS', 'RGB', 'HEB', 'BHEB', 'RHEB', 'EAGB', 'TPAGB']
        self.contam_header = '# %s %s \n' % (' '.join(self.regions),'Total')

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

        for i in range(len(self.sfr_files)):
            lines[-3] = ' '.join([self.sfr_files[i], extra])
            if object_mass is not None:
                lines[-6] = ' '.join(['%.4e' % object_mass, extra2])
            new_name = \
                os.path.split(self.galaxy_input)[1].replace('.dat',
                                                            '_%003i.dat' % i)
            new_out = os.path.join(self.outfile_loc, new_name)
            if dry_run is False:
                with open(new_out, 'w') as f:
                    f.write(''.join(lines))
                logger.info('wrote %s' % new_out)
            self.galaxy_inputs.append(new_out)

    def gather_results(self, mass_met=True, tpagb_lf=True, narratio_dict=None):
        '''gather results into strings or lists of strings for writing.'''
        result_dict = {}

        if tpagb_lf is True:
            # load mags
            if not hasattr(self, 'opt_mag'):
                self.load_trilegal_data()

            opt_hist, opt_bins = np.histogram(self.opt_mag, bins=self.opt_bins)

            ir_hist, ir_bins = np.histogram(self.ir_mag, bins=self.ir_bins)

            # scale the simulated LF to match the data LF
            opt_hist = np.array(opt_hist, dtype=float) * self.opt_norm
            ir_hist = np.array(ir_hist, dtype=float) * self.ir_norm

            result_dict['opt_lf_line'] = \
                '\n'.join([' '.join(['%g' % t for t in opt_hist]),
                           ' '.join(['%g' % t for t in opt_bins[1:]])])

            result_dict['ir_lf_line'] = \
                '\n'.join([' '.join(['%g' % t for t in ir_hist]),
                           ' '.join(['%g' % t for t in ir_bins[1:]])])

        if mass_met is True:
            self.sgal.all_stages('TPAGB')
            opt_inds = self.sgal.itpagb
            ir_inds = self.sgal.itpagb

            mag2 = self.sgal.mag2
            mag4 = self.sgal.data.get_col('F160W')

            opt_key = 'opt_mass_met_line'
            ir_key = 'ir_mass_met_line'
            for mag, inds, key in zip([mag2, mag4], [opt_inds, ir_inds],
                                      [opt_key, ir_key]):
                mass = self.sgal.data.get_col('m_ini')[inds]
                mag = mag[inds]
                mh = self.sgal.data.get_col('[M/H]')[inds]
                result_dict[key] = \
                    '\n'.join([' '.join(['%g' % t for t in mag]),
                               ' '.join(['%g' % t for t in mass]),
                               ' '.join(['%.3f' % t for t in mh])])


        # N agb/rgb ratio file
        narratio_dict = narratio_dict or {}
        if len(narratio_dict) > 0:
            self.narratio_header = '# target nopt_rgb nopt_agb nir_rgb nir_agb '
            self.narratio_header += 'opt_ar_ratio ir_ar_ratio opt_ar_ratio_err '
            self.narratio_header += 'ir_ar_ratio_err \n'

            narratio_fmt = '%(target)s %(nopt_rgb)i %(nopt_agb)i %(nir_rgb)i '
            narratio_fmt += '%(nir_agb)i %(opt_ar_ratio).3f %(ir_ar_ratio).3f '
            narratio_fmt += '%(opt_ar_ratio_err).3f  %(ir_ar_ratio_err).3f'

            opt_rgb = narratio_dict['opt_rgb']
            opt_agb = narratio_dict['opt_agb']
            ir_rgb = narratio_dict['ir_rgb']
            ir_agb = narratio_dict['ir_agb']

            nopt_rgb = float(len(opt_rgb))
            nopt_agb = float(len(opt_agb))
            nir_rgb = float(len(ir_rgb))
            nir_agb = float(len(ir_agb))
            out_dict = {'target': self.target,
                        'opt_ar_ratio': nopt_agb / nopt_rgb,
                        'ir_ar_ratio': nir_agb / nir_rgb,
                        'opt_ar_ratio_err':
                            galaxy_tests.count_uncert_ratio(nopt_agb, nopt_rgb),
                        'ir_ar_ratio_err':
                            galaxy_tests.count_uncert_ratio(nir_agb, nir_rgb),
                        'nopt_rgb': nopt_rgb,
                        'nopt_agb': nopt_agb,
                        'nir_rgb': nir_rgb,
                        'nir_agb': nir_agb}
            result_dict['narratio_line'] = narratio_fmt % out_dict

        return result_dict

    def do_normalization(self, filter1=None, trilegal_output=None,
                         hist_it_up=False, dry_run=False, debug=False):
        '''Do the normalization and save small part of outputs.'''
        if not hasattr(self, 'sgal') or self.mc is True:
            assert trilegal_output is not None, \
                'need sgal loaded or pass trilegal catalog file name'
            if filter1 is None:
                filter1 = self.filter1
            self.read_trilegal_catalog(trilegal_output, filter1=filter1)
            self.load_trilegal_data()

        # define regions
        opt_low = self.opt_offset
        opt_mid = self.opt_trgb + self.opt_trgb_err * self.ags.factor[0]
        opt_high = 10.

        ir_low = self.ir_offset
        ir_mid = self.ir_trgb + self.ir_trgb_err * self.ags.factor[1]
        ir_high = 10.

        # Recovered stars in simulated RGB region.
        sopt_rgb = self.sgal.stars_in_region(self.opt_mag, opt_low, opt_mid)
        sir_rgb = self.sgal.stars_in_region(self.ir_mag, ir_low, ir_mid)

        opt_mid = self.opt_trgb - self.opt_trgb_err * self.ags.factor[0]
        ir_mid = self.ir_trgb - self.ir_trgb_err * self.ags.factor[1]

        # Recovered stars in simulated AGB region.
        sopt_agb = self.sgal.stars_in_region(self.opt_mag, opt_mid, opt_high)
        sir_agb = self.sgal.stars_in_region(self.ir_mag, ir_mid, ir_high)

        # normalization
        self.opt_norm = self.nopt_rgb / float(len(sopt_rgb))
        self.ir_norm = self.nir_rgb / float(len(sir_rgb))

        logger.info('OPT Normalization: %f' % self.opt_norm)
        logger.info('IR Normalization: %f' % self.ir_norm)

        # random sample the data distribution
        rands = np.random.random(len(self.opt_mag))
        opt_ind, = np.nonzero(rands < self.opt_norm)
        rands = np.random.random(len(self.ir_mag))
        ir_ind, = np.nonzero(rands < self.ir_norm)

        # scaled rgb: norm + in rgb
        opt_rgb = list(set(opt_ind) & set(sopt_rgb))
        ir_rgb = list(set(ir_ind) & set(sir_rgb))

        # scaled agb
        opt_agb = list(set(opt_ind) & set(sopt_agb))
        ir_agb = list(set(ir_ind) & set(sir_agb))

        narratio_dict = {'opt_rgb': opt_rgb, 'opt_agb': opt_agb,
                        'ir_rgb': ir_rgb, 'ir_agb': ir_agb}
        result_dict = self.gather_results(narratio_dict=narratio_dict)
        if debug is True:
            return narratio_dict, (sopt_rgb, sopt_agb, sir_rgb, sir_agb), result_dict
        else:
            return (sopt_rgb, sopt_agb, sir_rgb, sir_agb), result_dict

    def vary_the_SFH(self, make_many_kw=None, hist_it_up=False):
        '''
        make the sfhs, make the galaxy inputs, run trilegal. For no trilegal
        runs, set dry_run True.
        '''

        make_many_kw = make_many_kw or {}
        if not 'mk_tri_sfh_kw' in make_many_kw.keys():
            make_many_kw['mk_tri_sfh_kw'] = {}

        dry_run = make_many_kw['mk_tri_sfh_kw'].get('dry_run', False)
        make_many_kw = dict({'nsfhs': self.nsfhs}.items() + \
                            make_many_kw.items())

        new_fmt = self.target + '_tri_%003i.sfr'
        sfr_outfilefmt = os.path.join(self.outfile_loc, new_fmt)

        make_many_kw['mk_tri_sfh_kw'] = \
            dict({'outfile_fmt': sfr_outfilefmt}.items() + \
                 make_many_kw['mk_tri_sfh_kw'].items())

        # what's up to the logger
        add_file_logger(self.outfile_loc)
        logger.debug(pprint.pformat(locals()))

        self.prepare_trilegal_sfr(make_many_kw=make_many_kw)

        self.prepare_galaxy_input(dry_run=dry_run)

        tname = os.path.join(self.outfile_loc,
                             'output_%s_%s' % (self.target, self.agb_mod))

        trilegal_output_fmt = tname + '_%003i.dat'

        result_dicts = []
        for galaxy_input in self.galaxy_inputs:
            num = int(galaxy_input.split('_')[-1].replace('.dat', ''))
            trilegal_output = trilegal_output_fmt % num
            rsp.TrilegalUtils.run_trilegal(self.cmd_input_file, galaxy_input,
                                           trilegal_output, rmfiles=False,
                                           dry_run=dry_run)

            norm_out, result_dict = \
                self.do_normalization(filter1=self.filter1,
                                      trilegal_output=trilegal_output,
                                      hist_it_up=hist_it_up,
                                      dry_run=dry_run)
            #self.binary_contamination(opt_agb, ir_agb)
            result_dict['contam_line'] = self.contamination_by_phases(*norm_out)
            result_dicts.append(result_dict)
            if dry_run is False:
                lastnum = int(num) - 1
                if os.path.isfile(trilegal_output_fmt % lastnum) is True:
                    os.remove(trilegal_output_fmt % lastnum)

        result = combine_list_of_dictionaries(result_dicts)
        return result

    def write_results(self, res_dict):
        '''writes out the results to self.fnames (see __init__)'''

        for fname in self.fnames:
            with open(fname, 'a') as fh:
                fshort = fname.split(self.target)[-1][1:].replace('.dat', '').replace(self.extra_str, '')
                if 'narratio' in fname:
                    fh.write(self.narratio_header)
                if 'contam' in fname:
                    fh.write(self.contam_header)

                [fh.write('%s \n' % l)for l in res_dict['%s_line' % fshort]]

        return

    def contamination_by_phases(self, sopt_rgb, sopt_agb, sir_rgb, sir_agb,
                                diag_plot=False):
        self.sgal.all_stages()
        indss = [self.sgal.__getattribute__('i%s' % r.lower()) for r in self.regions]
        line = ''
        contam_line = []
        if diag_plot is True:
            fig, (axs) = plt.subplots(ncols=2)
        for i, (rgb, agb, inds) in enumerate(zip([sopt_rgb, sir_rgb],
                                                 [sopt_agb, sir_agb],
                                                 [self.opt_color_cut,
                                                  self.ir_color_cut])):
            if i == 1:
                band = 'ir'
                mag = self.sgal.data.get_col('F160W')[inds]
            else:
                band = 'opt'
                mag = self.sgal.data.get_col('F814W')[inds]

            ncontam_rgb = [list(set(s) & set(inds) & set(rgb)) for s in indss]
            ncontam_agb = [list(set(s) & set(inds) & set(agb)) for s in indss]

            rheb_eagb_contam = len(ncontam_rgb[4]) + len(ncontam_rgb[5])
            frac_rheb_eagb = float(rheb_eagb_contam) / \
                float(np.sum([len(n) for n in ncontam_rgb]))

            heb_rgb_contam = len(ncontam_rgb[2])
            frac_heb_rgb_contam = float(heb_rgb_contam) / \
                float(np.sum([len(n) for n in ncontam_rgb]))

            mags = [mag[n] if len(n) > 0 else np.zeros(10) for n in ncontam_rgb]

            mms = np.concatenate(mags)
            ms, = np.nonzero(mms > 0)
            bins = np.linspace(np.min(mms[ms]), np.max(mms[ms]), 10)
            if diag_plot is True:
                [axs[i].hist(mags, bins=bins, alpha=0.5, stacked=True,
                             label=self.regions)]

            nrgb_cont = np.array([len(n) for n in ncontam_rgb], dtype=int)
            nagb_cont = np.array([len(n) for n in ncontam_agb], dtype=int)

            line += 'rgb %s %s %i \n' % (band, ' '.join(map(str, nrgb_cont)),
                                        np.sum(nrgb_cont))
            line += 'agb %s %s %i \n' % (band, ' '.join(map(str, nagb_cont)),
                                         np.sum(nagb_cont))

            line += '# rgb eagb contamination: %i \n' % rheb_eagb_contam
            line += '# frac of total in rgb region: %.3f \n' % frac_rheb_eagb
            line += '# rc contamination: %i \n' % heb_rgb_contam
            line += '# frac of total in rgb region: %.3f \n' % \
                    frac_heb_rgb_contam

            logger.info(line)

        contam_line.append(line)

        if diag_plot is True:
            axs[0].legend(numpoints=1, loc=0)
            axs[0].set_title(self.target)
            plt.savefig('contamination_%s.png' % self.target, dpi=150)

        return line

    def binary_contamination(self, sopt_agb, sir_agb):
        '''fraction of binaries in trgb (not implemented)'''
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


class Plotting(object):
    def __init__(self, vSFH):
        self.files = FileIO()
        bands = ['opt', 'ir']
        keys = ['trgb', 'trgb_err', 'offset']
        key = ['%s_%s' % (b, k) for b, k in itertools.product(bands, keys)]
        key += ['target', 'agb_mod', 'ags', 'fnames']
        [self.__setattr__(k, vSFH.__getattribute__(k)) for k in key]
        self.opt_lf_file, = [f for f in self.fnames if 'opt_lf' in f]
        self.ir_lf_file, = [f for f in self.fnames if 'ir_lf' in f]
        self.narratio_file, = [f for f in self.fnames if 'narratio' in f]

    def plot_lf_file(self, opt_lf_file, ir_lf_file, axs=None, plt_kw=None,
                     opt_limit=None, ir_limit=None):
        '''needs work, but: plot the lf files.'''
        # set up the plot
        plt_kw = plt_kw or {}
        plt_kw = dict({'linestyle': 'steps-mid', 'color': 'black',
                       'alpha': 0.2}.items() + plt_kw.items())
        label = '$%s$' % os.path.split(opt_lf_file)[1].split('_')[2].upper()
        plt_kw_lab = dict(plt_kw.items() + {'label': label}.items())
        if axs is None:
            fig, (axs) = plt.subplots(ncols=2, figsize=(12, 6))
            plt.subplots_adjust(right=0.95, left=0.05, wspace=0.1)

        # these have like 50 histograms each
        opt_hists, opt_binss = self.files.load_lf_file(self.opt_lf_file)
        ir_hists, ir_binss = self.files.load_lf_file(self.ir_lf_file)

        for i, (hists, binss, limit) in enumerate(zip([opt_hists, ir_hists],
                                                      [opt_binss, ir_binss],
                                                      [opt_limit, ir_limit])):

            for j, (hist, bins) in enumerate(zip(hists, binss)):
                if j != 0:
                    kw = plt_kw
                else:
                    kw = plt_kw_lab
                if limit is not None:
                    inds, = np.nonzero(bins <= limit)
                    axs[i].plot(bins[inds], hist[inds], **kw)
                else:
                    axs[i].plot(bins, hist, **kw)

        return axs

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
                                          trgb - trgb_err * self.ags.factor[i],
                                          10.)

            nrgb = np.sum(hist[irgb])
            nagb = np.sum(hist[iagb])
            ratio_data['%s_ar_ratio' % band] = nagb / nrgb
            ratio_data['%s_ar_ratio_err' % band] = \
                galaxy_tests.count_uncert_ratio(nagb, nrgb)
            ratio_data['n%s_rgb' % band] = nrgb
            ratio_data['n%s_agb'% band] = nagb
        return ratio_data

    def add_narratio_to_plot(self, ax, band, ratio_data):
        stext_kw = {'color': 'black', 'fontsize': 14, 'ha': 'center'}
        dtext_kw = {'color': 'darkred', 'fontsize': 14, 'ha': 'center'}
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

        textss = [[sagb_text, srgb_text, sratio_text],
                 [dagb_text, drgb_text, dratio_text]]
        kws = [stext_kw, dtext_kw]

        for kw, texts in zip(kws, textss):
            for xval, text in zip(xvals, texts):
                ax.text(xval, yval, text, **kw)
            yval -= .05  # stack the text
        return ax

    def plot_by_stage(self, ax1, ax2, add_stage_lfs='default', stage_lf_kw=None,
                      cols=None, trilegal_output=None, hist_it_up=False,
                      narratio=True):

        if add_stage_lfs == 'all':
            add_stage_lfs = ['RGB', 'HEB', 'RHEB',
                             'BHEB', 'EAGB', 'TPAGB']
        if add_stage_lfs == 'default':
            add_stage_lfs = ['RGB', 'EAGB', 'TPAGB']

        nstages = len(add_stage_lfs)
        stage_lf_kw = stage_lf_kw or {}
        stage_lf_kw = dict({'linestyle': 'steps', 'lw': 2}.items() +
                            stage_lf_kw.items())
        if cols is None:
            if nstages < 3:
                cmap = brewer2mpl.get_map('Paired', 'Qualitative', 3)
                cols = cmap.mpl_colors[0::2]
            else:
                cmap = brewer2mpl.get_map('Paired', 'Qualitative', nstages)
                cols = cmap.mpl_colors

        # load the trilegal catalog if it is given, if it is given,
        # no LF scaling... need to save this info better. Currently only
        # in log files.
        if trilegal_output is not None:
            self.files.read_trilegal_catalog(trilegal_output,
                                             filter1=get_filter1(self.target))
            self.files.load_trilegal_data()
            self.opt_norm = 1.
            self.ir_norm = 1.

        self.files.load_data_for_normalization(target=self.target, ags=self.ags)
        assert hasattr(self.files, 'opt_mag'), \
            'Need opt_mag or trilegal_output'

        for ax, mag, norm, sinds, bins in \
            zip([ax1, ax2],
                [self.files.sgal.data.get_col('F814W')-self.files.opt_moffset,
                 self.files.sgal.data.get_col('F160W')-self.files.ir_moffset],
                [self.opt_norm, self.ir_norm],
                [self.files.opt_color_cut, self.files.ir_color_cut],
                [self.files.opt_bins, self.files.ir_bins]):

            self.files.sgal.make_lf(mag, stages=add_stage_lfs, bins=bins,
                                    inds=sinds, hist_it_up=hist_it_up)
            for i in range(nstages):
                istage = add_stage_lfs[i].lower()
                try:
                    hist = self.files.sgal.__getattribute__('i%s_lfhist' %
                                                            istage)
                except AttributeError:
                    continue
                # combine all HeB stages into one for a cleaner plot.
                if add_stage_lfs[i].lower() == 'heb':
                    hist = \
                    np.sum([self.files.sgal.__getattribute__('i%s_lfhist' %
                                                             istage)
                            for j in range(nstages)
                            if 'heb' in istage], axis=0)
                elif 'heb' in istage:
                    continue

                bins = \
                self.files.sgal.__getattribute__('i%s_lfbins' %
                                                 istage)
                stage_lf_kw['color'] = cols[i]
                stage_lf_kw['label'] = '$%s$' % istage.upper()
                ax.plot(bins[:-1], hist*norm, **stage_lf_kw)

        sopt_hist, sopt_bins = self.files.sgal.make_lf(self.files.opt_mag,
                                                      bins=self.files.opt_bins,
                                                      hist_it_up=hist_it_up)
        sir_hist, sir_bins = self.files.sgal.make_lf(self.files.ir_mag,
                                                     bins=self.files.ir_bins,
                                                     hist_it_up=hist_it_up)

        sopt_hist = sopt_hist * self.opt_norm
        sir_hist = sir_hist * self.ir_norm
        stage_lf_kw['color'] = 'grey'

        lab = '_Total'
        if hasattr(self, 'agb_mod'):
            lab = self.agb_mod
        if lab != '':
            lab = '$%s$' % '\ '.join(lab.split('_')[1:])
        if narratio is False:
            stage_lf_kw['label'] = lab
        ax1.plot(sopt_bins[:-1], sopt_hist, **stage_lf_kw)
        ax2.plot(sir_bins[:-1], sir_hist, **stage_lf_kw)
        return ax1, ax2

    def compare_to_gal(self, hist_it_up=False, narratio=True, no_agb=False,
                       add_stage_lfs=None, extra_str='', trilegal_output=None,
                       plot_data=True, cols=None, stage_lf_kw=None, axs=None,
                       plt_kw=None, plot_models=True):
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
        # load plot limits:
        plims = load_plot_limits()
        # load ast_table for annotations
        ast_table = tables.read_completeness_table()

        # load galaxy data
        opt_gal, ir_gal = self.files.load_galaxies(hist_it_up=hist_it_up,
                                                   target=self.target,
                                                   ags=self.ags,
                                                   color_cut=True)

        if plot_models is True:
            # plot lfs from simulations (and initialize figure)
            plt_kw = plt_kw or {}
            (ax1, ax2) = \
                self.plot_lf_file(self.opt_lf_file, self.ir_lf_file,
                                  opt_limit=opt_gal.comp50mag2,
                                  ir_limit=ir_gal.comp50mag2, axs=axs,
                                  plt_kw=plt_kw)
        else:
            fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(12, 6),
                                           sharey=False)
            plt.subplots_adjust(right=0.95, left=0.1, wspace=0.1)

        if add_stage_lfs is not None and plot_models is False:
            (ax1, ax2) = self.plot_by_stage(add_stage_lfs, narratio=narratio,
                                            stage_lf_kw=stage_lf_kw, cols=cols,
                                            trilegal_output=trilegal_output)

        # plot galaxy data
        if plot_data is True:
            dplot_kw = \
                {'drawstyle': 'steps-mid', 'color': 'darkred', 'lw': 2,
                 'label': '$%s$' % self.target.upper().replace('-DEEP', '')}
            # HACK to mask low values
            #opt_gal.hist[opt_gal.hist < 0.1] = 0.1
            #ir_gal.hist[ir_gal.hist < 0.1] = 0.1
            opt_err = np.sqrt(opt_gal.hist)
            ir_err = np.sqrt(ir_gal.hist)
            ax1.errorbar(opt_gal.bins[1:], opt_gal.hist, yerr=opt_err,
                         **dplot_kw)
            ax2.errorbar(ir_gal.bins[1:], ir_gal.hist, yerr=ir_err, **dplot_kw)

        # initialize add numbers to the plot
        if narratio is True and plot_models is True:
            # count stars from the saved file
            ratio_data = rsp.fileIO.readfile(self.narratio_file,
                                             string_column=0)
            # get the number ratios for the annotations

            mean_ratio = {}
            for key in ratio_data.dtype.names:
                if key == 'target':
                    continue
                mean_ratio[key] = np.mean([ratio_data[i][key]
                                           for i in range(len(ratio_data))])

        for i, (ax, gal, trgb_err, band) in enumerate(zip([ax1, ax2],
                                                      [opt_gal, ir_gal],
                                                      [self.opt_trgb_err,
                                                      self.ir_trgb_err],
                                                      ['opt', 'ir'])):
            ax.set_yscale('log')
            ax.set_ylim(3, ax.get_ylim()[1])
            # set the max to be the brightest of the data .. not great.
            index = list(plims['target']).index(opt_gal.target)
            xmax = plims[index]['%s_cmdmin' % band]
            ax.set_xlim(xmax, gal.comp50mag2)

            yarr = np.linspace(*ax.get_ylim())
            # vertical lines around the trgb exclude region
            ax.fill_betweenx(yarr, gal.trgb - trgb_err * self.ags.factor[i],
                             gal.trgb + trgb_err * self.ags.factor[i],
                             color='black', alpha=0.2)

            ax.vlines(gal.trgb, *ax.get_ylim(), color='black',
                      linestyle='--')

            # % completeness limit
            ast_frac = rsp.fileIO.item_from_row(ast_table, 'target',
                                                self.target.upper(),
                                                '%s_filter2' % band)
            # line at dim limit for rgb normalization
            ax.vlines(ast_frac, *ax.get_ylim(), linestyle='--',
                      color='black')

            ax.fill_betweenx(yarr, ast_frac, ax.get_xlim()[1],
                             color='black', alpha=0.2)
            loc = 4
            if narratio is False:
                loc = 0
            ax.legend(loc=loc, frameon=False)
            ax.set_xlabel('$%s$' % gal.filter2, fontsize=20)

            if narratio is True:
                # need to load the data nrgb and nagb, calculate the ratio
                # and error.
                self.add_narratio_to_plot(ax, band, mean_ratio)

        ax1.set_ylabel('$\#$', fontsize=20)
        plt.tick_params(labelsize=16)
        outfile = '%s%s_lfs.png' % (self.opt_lf_file.split('opt_lf')[0][:-1],
                                   extra_str)
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
        target = '_'.join(os.path.split(opt_mass_met_file)[1].split('_')[0:4])
        fig.suptitle('$%s$' % target.replace('_', '\ '), fontsize=20)
        plt.savefig('%s_mass_met%s.png' % (target, extra_str), dpi=150)
        return grid


if __name__ == '__main__':
    #paolas_tests()
    ags = AncientGalaxies()
    galaxy_table = ags.write_trgb_table(exclude_region=[0.1, 0.2],
                                        mag_below_trgb='comp_frac')
    targets = ['ddo78', 'ddo71', 'hs117', 'kkh37', 'ngc2976', 'ngc404']
    #targets = ['ngc2976-deep', 'ngc404-deep']
    #targets = ['ddo71']
    mk_tri_sfh_kw = {'random_sfr': True, 'random_z': False}
    #call_vary_sfhs(targets, ['cmd_input_CAF09_S_NOV13.dat',
    #                                    'cmd_input_CAF09_S_NOV13eta0.dat',
    #                                    'cmd_input_CAF09_S_OCT13.dat'],
    #                          2, mk_tri_sfh_kw=mk_tri_sfh_kw, debug=False, dry_run=False)
