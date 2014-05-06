from star_formation_histories import StarFormationHistories
import os
from pop_synth.stellar_pops import normalize_simulation, rgb_agb_regions
import ResolvedStellarPops as rsp
import numpy as np
import matplotlib.pylab as plt

__all__ = ['setup_files', 'VarySFHs']

def setup_files(agb_mod, target, outfile_loc, extra_str=''):

    names = ['opt_lf', 'ir_lf', 'narratio', 'opt_mass_met', 'ir_mass_met',
             'contam']
    name_fmt = '%s_%s_%s%s.dat'

    fnames = [os.path.join(outfile_loc,
                           name_fmt % (agb_mod, target, f, extra_str))
              for f in names]

    return fnames


class VarySFHs(StarFormationHistories):
    '''
    run several variations of the age sfr z from MATCH SFH to produce
    simulated CMDs, LFs, and nagb/nrgb ratios.

    because the trilegal output can be so huge, there are many constraints
    needed to cull the model output down to what is necessary for the
    analysis.
    '''
    def __init__(self, **kwargs):
        '''
        galaxy_input is a template.
        '''
        # load SFH instance to make lots of trilegal runs

        default_kwargs = dict({'agb_mod': None,
                               'Av': 0.,
                               'cmd_input_file': None,
                               'dmod': 0.,
                               'extra_str': '',
                               'file_origin': None,
                               'filter1': None,
                               'galaxy_input': None,
                               'ir_bins': None,
                               'ir_color_min': -99.,
                               'ir_trgb': None,
                               'just_once': False,
                               'Mtrgb': None,
                               'nir_rgb': None,
                               'nopt_rgb': None,
                               'nsfhs': None,
                               'offsets': [-2., -1.5],
                               'opt_bins': None,
                               'opt_color_min': -99.,
                               'opt_trgb': None,
                               'outfile_loc': os.getcwd(),
                               'photsys': None,
                               'sfh_file': None,
                               'target': None,
                               'trgb_excludes': [.1, .2]}.items()
                               + kwargs.items())

        [self.__setattr__(k, v) for k, v in default_kwargs.items()]

        if self.just_once is False:
            StarFormationHistories.__init__(self, self.sfh_file, self.file_origin)

        if self.agb_mod is None:
            self.agb_mod = \
                os.path.split(self.cmd_input_file)[1].replace('.dat', '').lower()

        if None in [self.ir_trgb, self.opt_trgb]:
            assert self.Mtrgb is not None, 'need to supply trgb data'
            self.ir_trgb = rsp.astronomy_utils.Mag2mag(self.Mtrgb, 'F160W',
                                                  self.photsys, dmod=self.dmod,
                                                  Av=self.Av)
            self.opt_trgb = rsp.astronomy_utils.Mag2mag(self.Mtrgb, 'F814W',
                                                   self.photsys, dmod=self.dmod,
                                                   Av=self.Av)

        # setup the locations all the files to write and read from
        self.fnames =  setup_files(self.agb_mod, self.target, self.outfile_loc,
                                   extra_str=self.extra_str)

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

    def prepare_galaxy_input(self, object_mass=None, dry_run=False):
        '''
        write the galaxy input file from a previously written template.
        simply overwrites the filename line to link to the new sfr
        file.
        '''
        self.galaxy_inputs = []
        ext = '.' + self.galaxy_input.split('.')[-1]
        lines = open(self.galaxy_input).readlines()
        # line that links to sfr file.
        extra = ' '.join(lines[-3].split(' ')[1:])

        if object_mass is not None:
            extra2 = ' '.join(lines[-6].split()[1:])

        for i in range(len(self.sfr_files)):
            lines[-3] = ' '.join([self.sfr_files[i], extra])
            if object_mass is not None:
                lines[-6] = ' '.join(['%.4e' % object_mass, extra2]) + '\n'
            new_name = \
                os.path.split(self.galaxy_input)[1].replace(ext,
                                                            '_%003i' % i + ext )
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
                            rsp.math.utils.count_uncert_ratio(nopt_agb, nopt_rgb),
                        'ir_ar_ratio_err':
                            rsp.math.utils.count_uncert_ratio(nir_agb, nir_rgb),
                        'nopt_rgb': nopt_rgb,
                        'nopt_agb': nopt_agb,
                        'nir_rgb': nir_rgb,
                        'nir_agb': nir_agb}
            result_dict['narratio_line'] = narratio_fmt % out_dict

        return result_dict

    def load_trilegal_catalog(self, trilegal_output, shift_mags=False):
        '''read the trilegal cat mag1 and mag2 are optical.'''
        self.sgal = rsp.galaxies.SimGalaxy(trilegal_output,
                                           filter1=self.filter1,
                                           filter2='F814W')

        opt_mag = self.sgal.data.get_col('F814W')
        ir_mag = self.sgal.data.get_col('F160W')
        opt_mag1 = self.sgal.mag1
        ir_mag1 = self.sgal.data.get_col('F110W')

        opt_color_cut, = np.nonzero((opt_mag1 - opt_mag) > self.opt_color_min)
        ir_color_cut, = np.nonzero((ir_mag1 - ir_mag) > self.ir_color_min)


        if shift_mags is True:
            # Threshold for number of rgb stars in a bin.
            rgb_thresh = 5.
            rgb_bins = np.arange(10, 30, 0.01)
            self.sgal.all_stages('RGB')

            rgb_inds = np.intersect1d(self.sgal.irgb, opt_color_cut)
            rgb_hist, _ = np.histogram(opt_mag[rgb_inds], bins=rgb_bins)
            rgb_bin_edge = np.nonzero(rgb_hist > rgb_thresh)[0][0] - 1
            opt_offset = rgb_bins[rgb_bin_edge] - self.opt_trgb

            rgb_inds = np.intersect1d(self.sgal.irgb, ir_color_cut)
            rgb_hist, _ = np.histogram(ir_mag[rgb_inds], bins=rgb_bins)
            rgb_bin_edge = np.nonzero(rgb_hist > rgb_thresh)[0][0] - 1
            ir_offset = rgb_bins[rgb_bin_edge] - self.ir_trgb
        else:
            ir_offset = 0.
            opt_offset = 0.

        self.ir_moffset = ir_offset
        self.opt_moffset = opt_offset
        self.opt_color_cut = opt_color_cut
        self.ir_color_cut = ir_color_cut
        self.opt_mag = opt_mag[opt_color_cut] - opt_offset
        self.ir_mag = ir_mag[ir_color_cut] - ir_offset
        return

    def do_normalization(self, triout=None, debug=False):
        '''Do the normalization and save small part of outputs.'''
        if not hasattr(self, 'sgal') or self.nsfhs > 1:
            assert triout is not None, \
                'need sgal loaded or pass trilegal catalog file name'
            self.load_trilegal_catalog(triout)

        # select rgb and agb regions
        sopt_rgb, sir_rgb, sopt_agb, sir_agb = \
            rgb_agb_regions(self.sgal, self.offsets,
                            self.trgb_excludes, self.opt_trgb,
                            self.ir_trgb, self.opt_mag,
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

    def vary_the_SFH(self, random_sfr=True, random_z=False,
                     zdisp=True, dry_run=False, object_mass=None):
        '''
        make the sfhs, make the galaxy inputs, run trilegal. For no trilegal
        runs, set dry_run True.
        '''
        new_fmt = self.target + '_tri_%003i.sfr'
        outfile_fmt = os.path.join(self.outfile_loc, new_fmt)
        self.sfr_files = self.make_many_trilegal_sfhs(nsfhs=self.nsfhs,
                                                      outfile_fmt=outfile_fmt,
                                                      random_sfr=random_sfr,
                                                      random_z=random_z,
                                                      zdisp=zdisp,
                                                      dry_run=dry_run)

        self.prepare_galaxy_input(dry_run=dry_run, object_mass=object_mass)

        tname = os.path.join(self.outfile_loc,
                             'output_%s_%s' % (self.target, self.agb_mod))
        triout_fmt = tname + '_%003i.dat'

        result_dicts = []
        for i, galaxy_input in enumerate(self.galaxy_inputs):
            triout = triout_fmt % i

            rsp.trilegal.utils.run_trilegal(self.cmd_input_file, galaxy_input,
                                           triout, rmfiles=False,
                                           dry_run=dry_run)

            norm_out, result_dict = self.do_normalization(triout=triout)

            result_dict['contam_line'] = self.contamination_by_phases(*norm_out)
            result_dicts.append(result_dict)
            # remove the last trilegal output to save space
            if dry_run is False:
                lastnum = i - 1
                if os.path.isfile(triout_fmt % lastnum) is True:
                    os.remove(triout_fmt % lastnum)

        result = rsp.tools.helpers.combine_list_of_dictionaries(result_dicts)
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

    def load_lf_file(self, lf_file):
        with open(lf_file, 'r') as lff:
            lines = [l.strip() for l in lff.readlines()
                     if not l.startswith('#')]

        hists = [np.array(l.split(), dtype=float) for l in lines[0::2]]
        binss = [np.array(l.split(), dtype=float) for l in lines[1::2]]
        return hists, binss

def narratio_table(self):
    narratio_files = rsp.fileIO.get_files(self.outfile_dir, '*narratio*dat')
    stats.narratio_table(narratio_files)
    return

def chi2_stats(targets, cmd_inputs, outfile_dir='default', extra_str=''):
    chi2_files = stats.write_chi2_table(targets, cmd_inputs,
                                            outfile_loc=outfile_dir,
                                            extra_str=extra_str)
    chi2_dicts = stats.result2dict(chi2_files)
    stats.chi2plot(chi2_dicts, outfile_loc=outfile_dir)
    chi2_files = stats.write_chi2_table(targets, cmd_inputs,
                                            outfile_loc=outfile_dir,
                                            extra_str=extra_str,
                                            just_gauss=True)
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
