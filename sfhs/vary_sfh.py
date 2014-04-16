from star_formation_histories import StarFormationHistories
import os
from ..io import FileIO
from .. import tables
from ..TPAGBparams import phat_src
from ..math.utils import count_uncert_ratio
from ..utils import combine_list_of_dictionaries
from ..pop_synth.stellar_pops import normalize_simulation, rgb_agb_regions
import ResolvedStellarPops as rsp
import numpy as np
import matplotlib.pylab as plt


def setup_files(cmd_input_file, target, outfile_loc='default', extra_str='',
                mc=True, extra_directory='default'):

    agb_mod = None
    agb_model = cmd_input_file.replace('cmd_input_', '').lower()
    if extra_directory == 'default':
        agb_mod = agb_model.split('.')[0]

    if outfile_loc == 'default':
        outfile_loc = os.path.join(phat_src, 'Padova', 'with_agb', target)

        if agb_mod is not None:
            outfile_loc = os.path.join(outfile_loc, agb_mod)

        if mc is True:
            outfile_loc += '/mc'

        outfile_loc += '/'
        rsp.fileIO.ensure_dir(outfile_loc)

    names = ['opt_lf', 'ir_lf', 'narratio', 'opt_mass_met', 'ir_mass_met',
             'contam']
    name_fmt = '%s_%s_%s%s.dat'

    fnames = [os.path.join(outfile_loc,
                           name_fmt % (agb_mod, target, f, extra_str))
              for f in names]

    return outfile_loc, fnames, agb_mod


class VarySFHs(StarFormationHistories, FileIO):
    '''
    run several variations of the age sfr z from MATCH SFH to produce
    simulated CMDs, LFs, and nagb/nrgb ratios.

    can make plots, and save summary files.

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
            # e.g., fullpath/input_ddo68_stuff.dat
            gname = os.path.split(galaxy_input)[1]
            target = gname.split('_')[1].replace('.dat', '').lower()
        self.target = target
        self.filter1 = filter1
        self.galaxy_input = galaxy_input
        self.cmd_input_file = cmd_input_file
        self.nsfhs = nsfhs
        self.mc = False
        if self.nsfhs > 1:
            self.mc = True
        self.extra_str = extra_str

        # exclude regions and the number of data rgb and agb stars
        self.factor = [0., 0.]
        self.offests = 0.
        self.load_data_for_normalization()

        # 90% (or whatever) completeness magnitudes
        if len(self.offsets) == 0:
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
                print 'wrote %s' % new_out
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
                            count_uncert_ratio(nopt_agb, nopt_rgb),
                        'ir_ar_ratio_err':
                            count_uncert_ratio(nir_agb, nir_rgb),
                        'nopt_rgb': nopt_rgb,
                        'nopt_agb': nopt_agb,
                        'nir_rgb': nir_rgb,
                        'nir_agb': nir_agb}
            result_dict['narratio_line'] = narratio_fmt % out_dict

        return result_dict

    def do_normalization(self, filter1=None, trilegal_output=None,
                         hist_it_up=False, debug=False):
        '''Do the normalization and save small part of outputs.'''
        if not hasattr(self, 'sgal') or self.mc is True:
            assert trilegal_output is not None, \
                'need sgal loaded or pass trilegal catalog file name'
            if filter1 is None:
                filter1 = self.filter1
            self.read_trilegal_catalog(trilegal_output, filter1=filter1)
            self.load_trilegal_data()

        # select rgb and agb regions
        sopt_rgb, sir_rgb, sopt_agb, sir_agb = \
            rgb_agb_regions(self.sgal, self.opt_offset, self.opt_trgb,
                            self.opt_trgb_err, self.factor, self.ir_offset,
                            self.ir_trgb, self.ir_trgb_err, self.opt_mag,
                            self.ir_mag)

        # normalization
        self.opt_norm, self.ir_norm, opt_rgb, ir_rgb, opt_agb, ir_agb = \
            normalize_simulation(self.opt_mag, self.ir_mag, self.nopt_rgb,
                                 self.nir_rgb, sopt_rgb, sir_rgb, sopt_agb,
                                 sir_agb)

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

        new_fmt = self.target + '_tri_%003i.sfr'
        sfr_outfilefmt = os.path.join(self.outfile_loc, new_fmt)
        make_many_kw = dict({'nsfhs': self.nsfhs}.items() + \
                            make_many_kw.items())
        make_many_kw['mk_tri_sfh_kw'] = \
            dict({'outfile_fmt': sfr_outfilefmt}.items() + \
                 make_many_kw['mk_tri_sfh_kw'].items())

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
                                      hist_it_up=hist_it_up)
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
        indss = [self.sgal.__getattribute__('i%s' % r.lower())
                 for r in self.regions]
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

            print line

        contam_line.append(line)

        if diag_plot is True:
            axs[0].legend(numpoints=1, loc=0)
            axs[0].set_title(self.target)
            plt.savefig('contamination_%s.png' % self.target, dpi=150)

        return line


