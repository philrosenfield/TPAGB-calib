import logging
logging.basicConfig(filename='galaxy_tests.log',level=logging.DEBUG)
logger = logging.getLogger()
logger.info('start of run')
import ResolvedStellarPops as rsp
import ResolvedStellarPops.convertz as convertz

from mpl_toolkits.axes_grid1 import ImageGrid
from astroML.stats import binned_statistic_2d
import os
import numpy as np

import matplotlib.pyplot as plt
from TPAGBparams import table_src, snap_src
angst_data = rsp.angst_tables.AngstTables()
import galaxy_tests
import random

def number_of_stars(gal, exclude_region='default', mag_below_trgb=2.):
        '''
        Count the number of rgb and agb stars (in F814W or F160W)

        ARGS:
        exclude_region: mags above and below the trgb to exclude from
                        star counts [twice the trgb mag error]
        mag_below_trgb: mags below the trgb to call rgb stars [2]

        RETURNS:
        nrgb: int number of rgb stars
        nagb: int number of agb stars

        There are two reasons to use exclude_region.
        1) the TRGB mag is uncertain, and when counting agb stars,
           there should be no chance of scattering up rgb stars.
        2) TPAGB stars scatter below the TRGB due to the P in TP.
           Statistically, there are few TPAGB stars below the TRGB.

        TO DO:
        Test the number of AGB stars below the TRGB as a function of
        mag to test how many mTRGB_err exclude_region should use.

        Adapt for recent star formation:
        Add color information or verts to the stars in region call.
        Will need a robust contamination code.
        '''
        if exclude_region == 'default':
            if gal.filter2 == 'F814W':
                gal.trgb_err = angst_data.get_item(gal.target, 'mTRGB_err',
                                               extra_key='%s,%s' % (gal.filter1,
                                                                    gal.filter2))
            if gal.filter2 == 'F160W':
                gal.trgb_err = rsp.fileIO.item_from_row(angst_data.snap_tab3,
                                                    'target',
                                                    gal.target.upper(),
                                                    'mTRGB_F160W_err')

            exclude_region = gal.trgb_err * 2.

        # rgb will go from mag_below_trgb to trgb + exclude_region
        nrgb = len(gal.stars_in_region(gal.mag2, gal.trgb + mag_below_trgb,
                                       gal.trgb + exclude_region))

        # agb will go from trgb - exclude_region to very bright
        nagb = len(gal.stars_in_region(gal.mag2,
                                       gal.trgb - exclude_region, -99.))
        return nrgb, nagb, exclude_region


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
            tstring = targets
            targets = galaxy_tests.ancients()
        else:
            print 'need to write how to read other targets'

        fits_srcs = [snap_src + '/data/angst_no_trim', 'default']

        gal_dict = {}
        for i, band in enumerate(['opt', 'ir']):
            gals = rsp.Galaxies.galaxies([galaxy_tests.load_galaxy(t, band=band,
                                                      fits_src=fits_srcs[i])
                                          for t in targets])
            for gal in gals.galaxies:
                if 'ngc404-deep' == gal.target:
                    gal.target = 'ngc404'
                else:
                    gal.target = gal.target.lower()

                nrgb, nagb, ex_region = number_of_stars(gal,
                                                        mag_below_trgb=mag_below_trgb[i],
                                                        exclude_region=exclude_region)
                
                gdict = {'%s_trgb' % band: gal.trgb,
                         '%s_trgb_err' % band: gal.trgb_err,
                         'n%s_rgb' % band: nrgb,
                         'n%s_agb' % band: nagb}
                try:
                    gal_dict[gal.target].update(gdict)
                except KeyError:
                    gal_dict[gal.target] = gdict

        fmt = '%(target)s '
        fmt += '%(opt_trgb).2f %(opt_trgb_err).2f %(nopt_rgb)i %(nopt_agb)i '
        fmt += '%(ir_trgb).2f %(ir_trgb_err).2f %(nir_rgb)i %(nir_agb)i \n'
        outfile = os.path.join(table_src, '%s_galaxies.dat' % tstring)

        header = '# mags below trgb: optical=%.2f nir=%.2f \n' % (mag_below_trgb[0],
                                                                  mag_below_trgb[1])
        if exclude_region == 'default':
            header += '# excluded trgb region: %.2f * trgb mag_err \n' % ex_region
        header += '# target opt_trgb opt_trgb_err nopt_rgb nopt_agb '
        header += 'nir_trgb nir_trgb_err nir_rgb nir_agb \n'

        with open(outfile, 'w') as out:
            out.write(header)
            for k, v in gal_dict.items():
                gal_dict[k]['target'] = k
                out.write(fmt % gal_dict[k])
        print 'wrote %s' % outfile

    def read_trgb_table(self, table_name):
        dtype = [('target', '|S16'), ('opt_trgb', '<f8'),
                 ('opt_trgb_err', '<f8'), ('nopt_rgb', '<f8'),
                 ('nopt_agb', '<f8'), ('ir_trgb', '<f8'),
                 ('ir_trgb_err', '<f8'), ('nir_rgb', '<f8'),
                 ('nir_agb', '<f8')]
        data = np.genfromtxt(table_name, dtype=dtype)
        self.data = data.view(np.recarray)


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
    '''
    something
    '''
    def __init__(self, sfh_file, file_origin):
        self.base, self.name = os.path.split(sfh_file)
        self.data = parse_sfh_data(sfh_file, file_origin)

    def random_draw_within_uncertainty(self, attr, npoints=2e5):
        '''
        attr is the string name of the array that also has attr_errm and
        attr_errp (p and m are important due to the sign).

        If errm and errp are equal, just returns a randomly chosen point
        (of npoints) of a gaussian with mean attr and sigma = attr_errm

        If not, will stick to gaussians together at attr using sigma=attr_errm
        and sigm=attr_errp and returning a random value from there.

        If one of the err values is zero, will just use the other half of
        the gaussian.

        If they are both zero, well, just returns attr.
        '''
        #assert attr in ['sfr', 'mh'], 'Only set up for sfr and mh'
        if attr == 'mh':
            print 'Warning: this method was designed for sfr, not [M/H]'
        val_arr = self.data.__getattribute__(attr)
        errm_arr = self.data.__getattribute__('%s_errm' % attr)
        errp_arr = self.data.__getattribute__('%s_errp' % attr)
        rand_arr = np.array([])

        if attr == 'sfr':
            lowlim = 0
        else:
            lowlim = -999

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
                          zdisp=True, outfile='default'):
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
    This class was made to run several variations of the age sfr z
    zdisp table and save the resulting optical and nir LF, as well as
    the number of rgb and agb stars in the simulation.
    To save time, this need to be done first:
    table_file: Create a table of trgb and nrgb and nagb stars.
       Currently this is in AncientGalaxies and just takes all the
       stars within 1.5 or 2 mag below trgb. Will not be appropriate
       for galaxies with recent SF.
    
    the main method is vary_the_SFH, which sets up the files, calls
    StarFormationHistories and calls the other methods to write out the
    LFs and ratios.
    '''
    def __init__(self, galaxy_input, sfh_file, file_origin, target=None,
                 outfile_loc='default', table_file='default'):
        StarFormationHistories.__init__(self, sfh_file, file_origin)
        self.galaxy_input = galaxy_input

        if target is None:
            gname = os.path.split(self.galaxy_input)[1]
            target = gname.split('_')[1].replace('.dat', '').lower()
        self.target = target

        if outfile_loc == 'default':
            self.outfile_loc = os.path.join(snap_src, 'models', 'varysfh',
                                            self.target, 'mc/')

        rsp.fileIO.ensure_dir(self.outfile_loc)

        if table_file is 'default':
            table_file = os.path.join(table_src, 'ancients_galaxies.dat')
        # hmm better way?
        self.offsets = [2, 1.5]
        # else:
        # maybe write one?
        self.ags = AncientGalaxies()
        self.ags.read_trgb_table(table_file)

        # load stars in data RGB and AGB region as well as TRGBs.
        # this is just a call to self.ags.data
        self.load_data_for_normalization()

    def prepare_trilegal_sfr(self, make_many_kw=None):
        '''call make_many_trilegal_sfhs'''
        make_many_kw = make_many_kw or {}
        self.sfr_files = self.make_many_trilegal_sfhs(**make_many_kw)

    def prepare_galaxy_input(self, object_mass=None):
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
            
            with open(new_out, 'w') as f:
                f.write(''.join(lines))
            print 'wrote %s' % new_out
            self.galaxy_inputs.append(new_out)

    def prepare_outfiles(self):
        '''
        prepare outfiles for vary_the_sfh.
        NOTE: all attrs saved are file objects.
        '''
        # LF file
        opt_fname = os.path.join(self.outfile_loc,
                                 '%s_opt_lf.dat' % self.target)
        ir_fname = os.path.join(self.outfile_loc,
                                '%s_ir_lf.dat' % self.target)
        self.opt_lf_file = open(opt_fname, 'a')
        self.ir_lf_file = open(ir_fname, 'a')

        # N agb/rgb ratio file
        fmt =  '%(target)s %(nopt_rgb)i %(nopt_agb)i %(nir_rgb)i '
        fmt += '%(nir_agb)i %(opt_ar_ratio).3f %(ir_ar_ratio).3f '
        fmt += '%(opt_ar_ratio_err).3f  %(ir_ar_ratio_err).3f \n'
        
        header = '# target nopt_rgb nopt_agb nir_rgb nir_agb opt_ar_ratio '
        header += 'ir_ar_ratio opt_ar_ratio_err ir_ar_ratio_err \n'
        narratio_fname = os.path.join(self.outfile_loc,
                                      '%s_narratio.dat' % self.target)
        
        self.narratio_fmt = fmt
        self.narratio_file = open(narratio_fname, 'a')
        self.narratio_file.write(header)
        
        # mass_met_file
        omass_met_fname = os.path.join(self.outfile_loc,
                                      '%s_opt_mass_met.dat' % self.target)
        imass_met_fname = os.path.join(self.outfile_loc,
                                      '%s_ir_mass_met.dat' % self.target)

        self.opt_mass_met_file = open(omass_met_fname, 'a')
        self.ir_mass_met_file = open(imass_met_fname, 'a')
        
        # how to name the trilegal sfr files
        new_fmt = self.name + '_tri_%003i.sfr'
        self.sfr_outfilefmt = os.path.join(self.outfile_loc, new_fmt)
        return opt_fname, ir_fname, narratio_fname
    
    def close_files(self):
        '''close all open files'''
        [f.close() for f in self.__dict__.values() if type(f) is file]

    def write_LF(self, opt_norm, ir_norm, opt_file, ir_file):
        '''
        write the LF to files.
        ARGS:
        opt_file, ir_file: write out files, can be file objects.
        opt_rec and ir_rec are indices of sgal f814w f160w.

        NOTE: must have attr sgal loaded.
        Formatting is a repeat of first line is histogram, next line is bins[1:]
        '''
        # load mags, set here in case incase I want to pass this in the future.
        opt_mag = self.sgal.data.get_col('%s_cor' % self.ast_objs[0].filter2)
        ir_mag = self.sgal.data.get_col('%s_cor' % self.ast_objs[1].filter2)

        # hist it up!
        opt_hist, opt_bins = galaxy_tests.hist_it_up(opt_mag, threash=5)
        ir_hist, ir_bins = galaxy_tests.hist_it_up(ir_mag, threash=5)

        # scale the simulated LF to match the data LF
        opt_hist = np.array(opt_hist, dtype=float) * opt_norm
        ir_hist = np.array(ir_hist, dtype=float) * ir_norm

        np.savetxt(opt_file, (opt_hist, opt_bins[1:]), fmt='%g')
        np.savetxt(ir_file, (ir_hist, ir_bins[1:]), fmt='%g')

        print 'wrote %s, %s' % (opt_file.name, ir_file.name)
        return

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

        # recovered tpagb stars
        opt_inds = np.intersect1d(self.sgal.itpagb, self.sgal.rec)
        ir_inds = np.intersect1d(self.sgal.itpagb, self.sgal.ir_rec)

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
        fake_files = galaxy_tests.get_fake_files(self.target.upper())
        ast_objs = [rsp.Galaxies.artificial_star_tests(f)
                    for f in fake_files]
        self.ast_objs = ast_objs

    def read_trilegal_catalog(self, trilegal_output, filter1='F606W'):
        '''
        reads the trilegal cat and does ast corrections.
        '''
        sgal = rsp.Galaxies.simgalaxy(trilegal_output, filter1=filter1,
                                           filter2='F814W')
        if not hasattr(self, 'ast_objs'):
            # should be loaded outside this method if mc is running.
            self.load_asts()

        # ast correction
        rsp.Galaxies.ast_correct_trilegal_sim(sgal, asts_obj=self.ast_objs)
        sgal.load_ast_corrections()
        
        self.sgal = sgal

    def do_normalization(self, mc=True, filter1='F606W', trilegal_output=None):
        '''Do the normalization and save small part of outputs.'''
        if not hasattr(self, 'sgal'):
            assert trilegal_output is not None, \
                'need sgal loaded or pass trilegal catalog file name'
            self.read_trilegal_catalog(trilegal_output, filter1=filter1)

        opt_mag = self.sgal.data.get_col('%s_cor' % self.ast_objs[0].filter2)
        ir_mag = self.sgal.data.get_col('%s_cor' % self.ast_objs[1].filter2)

        # Recovered stars in simulated RGB region.
        sopt_rgb = self.sgal.stars_in_region(opt_mag,
                                             self.opt_trgb + self.offsets[0],
                                             self.opt_trgb + self.opt_trgb_err * 2)
        sir_rgb = self.sgal.stars_in_region(ir_mag,
                                            self.ir_trgb + self.offsets[1],
                                            self.ir_trgb + self.ir_trgb_err * 2)

        # Recovered stars in simulated AGB region.
        sopt_agb = self.sgal.stars_in_region(opt_mag, self.opt_trgb - self.opt_trgb_err * 2, 10.)
        sir_agb = self.sgal.stars_in_region(ir_mag, self.ir_trgb - self.ir_trgb_err * 2, 10.)

        # normalization
        opt_norm = self.nopt_rgb / float(len(sopt_rgb))
        ir_norm = self.nir_rgb / float(len(sir_rgb))

        #print opt_norm, ir_norm

        # random sample the data distribution
        rands = np.random.random(len(opt_mag))
        opt_ind, = np.nonzero(rands < opt_norm)
        rands = np.random.random(len(ir_mag))
        ir_ind, = np.nonzero(rands < ir_norm)

        # scaled rgb: ast + norm + in rgb
        opt_rgb = list(set(opt_ind) & set(sopt_rgb) & set(self.sgal.rec))
        ir_rgb = list(set(ir_ind) & set(sir_rgb) & set(self.sgal.ir_rec))

        #print len(opt_rgb), self.nopt_rgb
        #print len(ir_rgb), self.nir_rgb

        # scaled agb
        opt_agb = list(set(opt_ind) & set(sopt_agb) & set(self.sgal.rec))
        ir_agb = list(set(ir_ind) & set(sir_agb) & set(self.sgal.ir_rec))

        #save LF in both filters
        if mc is True:
            # only need to save the LF and the ratios.
            self.write_mass_met_file()
            self.write_LF(opt_norm, ir_norm, self.opt_lf_file, self.ir_lf_file)
            self.write_ratio(opt_rgb, opt_agb, ir_rgb, ir_agb,
                             self.narratio_file)
        return opt_norm, ir_norm, opt_ind, ir_ind

    def vary_the_SFH(self, cmd_input_file, make_many_kw=None, dry_run=False,
                     mc=True, diag_plots=False):
        '''
        make the sfhs, make the galaxy inputs, run trilegal. For no trilegal
        runs, set dry_run True.
        '''
        self.prepare_outfiles()
        
        make_many_kw = make_many_kw or {}
        if not 'mk_tri_sfh_kw' in make_many_kw.keys():
            make_many_kw['mk_tri_sfh_kw'] = {}

        make_many_kw['mk_tri_sfh_kw'] = dict({'outfile_fmt': self.sfr_outfilefmt}.items() +
                                              make_many_kw['mk_tri_sfh_kw'].items())
        self.prepare_trilegal_sfr(make_many_kw=make_many_kw)

        self.prepare_galaxy_input()

        agb_model = cmd_input_file.replace('cmd_input_', '').lower()
        tname = 'output_%s_%s' % (self.target, agb_model)
        trilegal_output = os.path.join(self.outfile_loc, tname)

        self.load_asts()
        filter1 = self.ast_objs[0].filter1

        for galaxy_input in self.galaxy_inputs:
            rsp.TrilegalUtils.run_trilegal(cmd_input_file, galaxy_input,
                                           trilegal_output, rmfiles=False,
                                           dry_run=dry_run)
            self.do_normalization(trilegal_output=trilegal_output,
                                  mc=mc, filter1=filter1)
        self.close_files()
        if diag_plots is True:
            outfile_fmt = os.path.join(self.outfile_loc, '%s_random' % self.target)
            outfile_fmt += '_%s.png'
            [self.plot_random_arrays(attr_str, from_files=True,
                                     outfile=outfile_fmt % attr_str)
             for attr_str in ['sfr', 'mh']]
            
            self.compare_to_gal()
            self.plot_mass_met_table(self.opt_mass_met_file.name,
                                     self.ir_mass_met_file.name)
    
    def plot_lf_file(self, opt_lf_file, ir_lf_file, axs=None):
        '''needs work, but: plot the lf files.'''
        if axs is None:
            fig, (axs) = plt.subplots(ncols=2, figsize=(8, 8))

        for i, lf_file in enumerate([opt_lf_file, ir_lf_file]):

            with open(lf_file, 'r') as lff:
                lines = [l.strip() for l in lff.readlines()
                         if not l.startswith('#')]

            hists = [np.array(l.split(), dtype=float) for l in lines[0::2]]
            binss = [np.array(l.split(), dtype=float) for l in lines[1::2]]

            for hist, bins in zip(hists, binss):
                if len(hist) != len(bins):
                    continue
                axs[i].plot(bins, hist, linestyle='steps', color='k', alpha=0.2)
        return axs

    def compare_to_gal(self):
        '''
        work in progress...
        '''
        ir_gal = galaxy_tests.load_galaxy(self.target, band='ir')
        ir_hist, ir_bins =  galaxy_tests.hist_it_up(ir_gal.mag2, threash=5)

        fits_src = snap_src + '/data/angst_no_trim'
        opt_gal = galaxy_tests.load_galaxy(self.target, band='opt',
                                           fits_src=fits_src)
        opt_hist, opt_bins =  galaxy_tests.hist_it_up(opt_gal.mag2, threash=5)

        ax1, ax2 = self.plot_lf_file(self.opt_lf_file.name,
                                     self.ir_lf_file.name)

        ax2.plot(ir_bins[1:], ir_hist, linestyle='steps', color='navy', lw=2)
        ax1.plot(opt_bins[1:], opt_hist, linestyle='steps', color='navy', lw=2,
                 label='$%s$' % self.target.upper())

        for ax, gal, offset in zip([ax1, ax2], [opt_gal, ir_gal], self.offsets):
            ax.set_yscale('log')
            ax.set_xlim(ax.get_xlim()[0], gal.comp50mag2)
            ax.vlines(gal.trgb, *ax1.get_ylim(), lw=2, color='darkred')
            ax.vlines(gal.trgb + offset, *ax.get_ylim(),
                      color='darkred')

            ax.set_xlabel('$%s$' % gal.filter2, fontsize=20)
            ax.set_ylabel('$\#$', fontsize=20)

        ax1.legend(loc=1, frameon=False)
        outfile = os.path.join(self.outfile_loc, '%s_lfs.png' % self.target)
        plt.savefig(outfile, dpi=300)
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
            for i, mass_met in enumerate([opt_mass_met_file, ir_mass_met_file]):
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
                
                N, xedges, yedges = binned_statistic_2d(mag, mass, mh, types[j], bins=50)
                im = grid[k].imshow(N.T, origin='lower',
                               extent=[xedges[0], xedges[-1], yedges[0], yedges[-1]],
                               aspect='auto', interpolation='nearest', cmap=cmaps[j])
                grid[k].cax.colorbar(im)
                #grid[i].cax.set_label('$[M/H]$')

        grid.axes_all[0].set_ylabel('${\\rm Mass}\ (M_\odot)$', fontsize=20)
        grid.axes_all[2].set_ylabel('${\\rm Mass}\ (M_\odot)$', fontsize=20)
        grid.axes_all[2].set_xlabel('$F814W$', fontsize=20)
        grid.axes_all[3].set_xlabel('$F160W$', fontsize=20)
        fig.suptitle('$DDO71$', fontsize=20)
        plt.savefig(os.path.join(self.outfile_loc, '%s_mass_met.png' % self.target), dpi=300)
        return grid

    def plot_ratio(self):
        from astroML.plotting import hist
        ratio = rsp.fileIO.readfile(self.narratio_file.name)
        hist(ratio['opt_ar_ratio'], bins='blocks')
        plt.hist(ratio['ir_ar_ratio'])
        pass

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
    # find match sfh files
    sfhs = rsp.fileIO.get_files(sfr_dir, search_str)

    # load SFHs
    SFHs = [StarFormationHistories(s, fileorigin) for s in sfhs]

    # write trilegal sfh
    tri_sfhs = [S.make_trilegal_sfh(**make_trilegal_sfh_kw) for S in SFHs]
    return tri_sfhs


def vary_sfhs_of_one_galaxy(galaxy_name, cmd_input_file, mk_tri_sfh_kw=None,
                            match_sfh_file='default', vary_sfh_kw=None,
                            match_fileorigin='match-grid', make_many_kw=None,
                            galaxy_input_file='default'):
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

    vary_sfh_kw = dict({'diag_plots': True,
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
                    target=galaxy_name.lower())

    # don't keep appending to old files.
    vSFH.remove_files()

    # vary the SFH
    vSFH.vary_the_SFH(cmd_input_file, **vary_sfh_kw)

    return vSFH
