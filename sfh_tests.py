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

    def number_of_stars(self, gal, offset, trgb_offset=0):
        '''
        more complex? Add color information or verts to the
        stars in region call. Also write a robust contamination code.
        '''
        nrgb = len(gal.stars_in_region(gal.mag2, gal.trgb + offset,
                                                gal.trgb + trgb_offset))
        nagb = len(gal.stars_in_region(gal.mag2, gal.trgb, 99.))
        return nrgb, nagb
        
    def write_trgb_table(self, targets='ancients', offsets=[2, 1.5],
                         trgb_offset=0):
        '''
        writes the trgb and number of stars within offsets (opt, ir) mags of
        trgb
        '''
        if targets == 'ancients':
            tstring = targets
            targets = galaxy_tests.ancients()
        else:
            print 'need to write how to read other targets'

        fits_srcs = [snap_src + '/data/angst_no_trim',
                     'default']
        gal_dict = {}
        for i, band in enumerate(['opt', 'ir']):
            print band
            gals = rsp.Galaxies.galaxies([galaxy_tests.load_galaxy(t, band=band,
                                                      fits_src=fits_srcs[i])
                                          for t in targets])
            for gal in gals.galaxies:
                if 'ngc404-deep' == gal.target:
                    gal.target = 'ngc404'
                else:
                    gal.target = gal.target.lower()

                nrgb, nagb = self.number_of_stars(gal, offsets[i])
                
                gdict = {'%s_trgb' % band: gal.trgb,
                         'n%s_rgb' % band: nrgb,
                         'n%s_agb' % band: nagb}
                try:
                    gal_dict[gal.target].update(gdict)
                except KeyError:
                    gal_dict[gal.target] = gdict

        fmt = '%(target)s %(opt_trgb).2f %(nopt_rgb)i %(nopt_agb)i %(ir_trgb).2f %(nir_rgb)i %(nir_agb)i \n'
        outfile = os.path.join(table_src, '%s_galaxies_%s_mag_below.dat' %
                              (tstring, '_'.join(map(str, offsets))))
        if trgb_offset > 0:
            outfile = outfile.replace('_below.dat', '_below_trgb_%s.dat' % str(trgb_offset))
            
        with open(outfile, 'w') as out:
            out.write('# target opt_trgb nopt_rgb nopt_agb ir_trgb nir_rgb nir_agb \n')
            for k, v in gal_dict.items():
                gal_dict[k]['target'] = k
                out.write(fmt % gal_dict[k])
        print 'wrote %s' % outfile

    def read_trgb_table(self, table_name):
        dtype = [('target', '|S16'), ('opt_trgb', '<f8'), ('nopt_rgb', '<f8'),
                 ('nopt_agb', '<f8'), ('ir_trgb', '<f8'), ('nir_rgb', '<f8'),
                 ('nir_agb', '<f8')]
        data = np.genfromtxt(table_name, dtype=dtype)
        self.data = data.view(np.recarray)


def parse_sfh_data(filename, file_origin):
    '''
    supply a match_sfh (.zs.sfh) file to parse, or np.recarray containing
        dtype = [('lagei', '<f8'),
                ('lagef', '<f8'),
                ('sfr', '<f8'),
                ('sfr_errp', '<f8'),
                ('sfr_errm', '<f8'),
                ('mh', '<f8'),
                ('mh_errp', '<f8'),
                ('mh_errm', '<f8'),
                ('mh_disp', '<f8')]
    '''
    if 'match' in file_origin.lower():
        data = rsp.match_utils.read_binned_sfh(filename)
    else:
        print 'please add a new data reader'
    
    if 'grid' in file_origin.lower():
        # at least have a uniform error that is 10% the smallest sfr.
        data.sfr_errp = np.min(data.sfr[data.sfr > 0]) * 0.2
        # now take 10% as nominal error in each sfr bin.
        data.sfr_errp[data.sfr > 0] = data.sfr[data.sfr > 0] * 0.2
        data.sfr_errm[data.sfr > 0] = data.sfr[data.sfr > 0] * 0.2
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
                                 self.name.replace('.zc.sfh', '.tri.dat'))

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
            val_arrs = np.array([10**(val_arr/.2)
                                 for val_arr in val_arrs])
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
                                       self.name.replace('.zc.sfh',
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
    This class was made to run several variations of the age sfr z zdisp
    table and save the resulting optical and nir LF, as well as the number of
    rgb and agb stars in the simulation.
    To save time, this need to be done first:
    table_file: Create a table of trgb and nrgb and nagb stars. Currently
       this is in AncientGalaxies and just takes all the stars within 1.5 or
       2 mag below trgb. Will not be appropriate for galaxies with recent SF.
    
    the main method is vary_the_SFH, which sets up the files, calls
    StarFormationHistories and calls the other methods to write out the LFs and
    ratios.
    '''
    def __init__(self, galaxy_input, sfh_file, file_origin, outfile_loc='default',
                 table_file='default', target=None):
        StarFormationHistories.__init__(self, sfh_file, file_origin)
        self.galaxy_input = galaxy_input

        if target is None:
            gname = os.path.split(self.galaxy_input)[1]
            target = gname.split('_')[1].replace('.dat', '').lower()
        self.target = target

        if outfile_loc == 'default':
            self.outfile_loc = os.path.join(snap_src, 'data', 'sfh_parsec',
                                            self.target, 'mc/')

        rsp.fileIO.ensure_dir(self.outfile_loc)

        if table_file is 'default':
            table_file = os.path.join(table_src,
                                      'ancients_galaxies_2_1.5_mag_below.dat')
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
        '''
        a call to make_many_trilegal_sfhs
        '''
        make_many_kw = make_many_kw or {}
        self.sfr_files = self.make_many_trilegal_sfhs(**make_many_kw)

    def prepare_galaxy_input(self, object_mass=None):
        '''
        writes the galaxy input file from a previously written template.
        simply overwrites the filename line to link to the new sfr file.
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
        note, all attrs saved are file objects. Probably should have something
        to not keep adding to the same files on different runs. 
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
        [f.close() for f in self.__dict__.values() if type(f) is file]

    def write_LF(self, opt_rec, ir_rec, opt_file, ir_file):
        '''
        writes out the LF to files. opt_file, ir_file can already be
        file objects.
        must have attr sgal loaded.
        opt_rec and ir_rec are indices of sgal f814w f160w.
        first line is histogram
        next line is bins[1:]
        etc.
        '''
        opt_mag = self.sgal.ast_mag2[opt_rec]
        ir_mag = self.sgal.ast_mag4[ir_rec]

        opt_hist, opt_bins = galaxy_tests.hist_it_up(opt_mag, threash=5)
        ir_hist, ir_bins = galaxy_tests.hist_it_up(ir_mag, threash=5)

        np.savetxt(opt_file, (opt_hist, opt_bins[1:]), fmt='%g')
        np.savetxt(ir_file, (ir_hist, ir_bins[1:]), fmt='%g')

        print 'wrote %s, %s' % (opt_file.name, ir_file.name)
        return

    def write_ratio(self, opt_rgb, opt_agb, ir_rgb, ir_agb, out_file):
        '''
        write the numbers of stars, the ratio, and the poisson uncertainty in
        the ratio.
        see prepare_files for fmt information
        out_file is an opened file object.
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

    def write_scaled_files(self):
        pass

    def write_mass_met_file(self):
        self.sgal.all_stages('TPAGB')
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
        '''
        load the numbers of data rgb and agb stars. The file was already read
        '''
        # need to get normalization...
        # for ancients, this is just a mag cut, no color cut.
        self.nopt_rgb = rsp.fileIO.item_from_row(self.ags.data, 'target',
                                                 self.target, 'nopt_rgb')
        self.nir_rgb = rsp.fileIO.item_from_row(self.ags.data, 'target',
                                                self.target, 'nir_rgb')
        self.opt_trgb = rsp.fileIO.item_from_row(self.ags.data, 'target',
                                                self.target, 'opt_trgb')
        self.ir_trgb = rsp.fileIO.item_from_row(self.ags.data, 'target',
                                                self.target, 'ir_trgb')
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
        '''
        Do the normalization and save small part of outputs.
        '''
        if not hasattr(self, 'sgal'):
            assert trilegal_output is not None, \
                'need sgal loaded or pass trilegal catalog file name'
            self.read_trilegal_catalog(trilegal_output, filter1='F606W')

        tget = self.target.upper()

        min_opt = angst_data.__getattribute__(tget)['F814W']['50_completeness']
        min_ir = angst_data.get_snap_50compmag(tget, 'F160W')

        imag2, = np.nonzero(self.sgal.ast_mag2 < min_opt)
        imag4, = np.nonzero(self.sgal.ast_mag4 < min_ir)

        mag2 = self.sgal.ast_mag2
        mag4 = self.sgal.ast_mag4

        # Recovered stars in simulated RGB region.
        sopt_rgb = self.sgal.stars_in_region(mag2,
                                             self.opt_trgb + self.offsets[0],
                                             self.opt_trgb)
        sir_rgb = self.sgal.stars_in_region(mag4,
                                            self.ir_trgb + self.offsets[1],
                                            self.ir_trgb)

        # Recovered stars in simulated AGB region.
        sopt_agb = self.sgal.stars_in_region(mag2, self.opt_trgb, 10.)
        sir_agb = self.sgal.stars_in_region(mag4, self.ir_trgb, 10.)

        # normalization
        opt_norm = self.nopt_rgb / float(len(sopt_rgb))
        sim_norm = self.nir_rgb / float(len(sir_rgb))

        print opt_norm, sim_norm

        # random sample the data distribution
        rands = np.random.random(len(self.sgal.ast_mag2))
        opt_ind, = np.nonzero(rands < opt_norm)
        ir_ind, = np.nonzero(rands < sim_norm)

        # scaled: ast + norm
        opt_rec = list(set(opt_ind) & set(self.sgal.rec) & set(imag2))
        ir_rec = list(set(ir_ind) & set(self.sgal.ir_rec) & set(imag4))

        # scaled rgb: ast + norm + in rgb
        opt_rgb = list(set(opt_ind) & set(sopt_rgb) & set(self.sgal.rec))
        ir_rgb = list(set(ir_ind) & set(sir_rgb) & set(self.sgal.ir_rec))

        print len(opt_rgb), self.nopt_rgb
        print len(ir_rgb), self.nir_rgb

        # scaled agb
        opt_agb = list(set(opt_ind) & set(sopt_agb) & set(self.sgal.rec))
        ir_agb = list(set(ir_ind) & set(sir_agb) & set(self.sgal.ir_rec))

        #save LF in both filters
        if mc:
            # only need to save the LF and the ratios.
            self.write_mass_met_file()
            self.write_LF(opt_rec, ir_rec, self.opt_lf_file, self.ir_lf_file)
            self.write_ratio(opt_rgb, opt_agb, ir_rgb, ir_agb,
                             self.narratio_file)
        else:
            # save scaled files to make full LFs and have choice in binning
            self.write_scaled_files(opt_rec, ir_rec)

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

        for galaxy_input in self.galaxy_inputs:
            rsp.TrilegalUtils.run_trilegal(cmd_input_file, galaxy_input,
                                           trilegal_output, rmfiles=False,
                                           dry_run=dry_run)
            self.do_normalization(trilegal_output=trilegal_output,
                                  mc=mc, filter1='F606W')
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
    
    def plot_lf_file(self, opt_lf_file, ir_lf_file):
        '''
        needs work, but will plot the lf files.
        '''
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

