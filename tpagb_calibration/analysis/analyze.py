"""
A dump of snippets from vary_sfh that shouldn't have been in there!
"""
import argparse
import logging
import numpy as np
import os
import sys
import time

import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pylab as plt
import ResolvedStellarPops as rsp

from IPython import parallel
from ..pop_synth.stellar_pops import normalize_simulation, rgb_agb_regions
from ..plotting.plotting import model_cmd_withasts
from star_formation_histories import StarFormationHistories

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)
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
def initialize_inputs():
    return {'ast_corr': True,
            'cmd_input_file': 'cmd_input_parsecCAF09_V1.2S_M36_S12D2.dat',
            'extra_str': '',
            'file_origin': None,
            'filter1': None,
            'filter2': None,
            'galaxy_input': None,
            'bins': None,
            'trgb': None,
            'Mtrgb': None,
            'offset': -2.,
            'outfile_loc': os.getcwd(),
            'photsys': None,
            'sfh_file': None,
            'target': None,
            'trgb_exclude': .1}


            """
            if do_norm:
                for r in res:
                    logging.info(r.result)
                    filter2 = self.filter2
                    filter1 = self.filter1
                    if self.ast_corr:
                        filter2 = '{}_cor'.format(self.filter2)
                        filter1 = '{}_cor'.format(self.filter1)
                    fdict = write_results(r.result, self.agb_mod, self.target,
                                         self.outfile_loc, filter2, filter1,
                                         extra_str=self.extra_str)
                    r.result.update(fdict)

            # to eliminate clutter
            if cleanup:
                for i in iset:
                    if i != self.nsfhs - 1:
                        if os.path.isfile(triout_fmt % i):
                            os.remove(triout_fmt % i)
            
            # write the new "input file"
            for i in range(len(res)):
                out_obj.add_params(res[i].result)
                out_obj.write_params(outparam % i, loud=True)
            """
        #os.system('ipcluster stop')

    def run(self, do_norm=True, dry_run=False, is_parallel=False, do_norm_kw={}):
        """call run_once and write results if normalization happened"""
        do_norm_kw = dict(self.__dict__.items() + do_norm_kw.items())
        try:
            del do_norm_kw['sgal']
        except:
            pass
        run_once(do_norm=do_norm, dry_run=dry_run, do_norm_kw=do_norm_kw)
        """
        self.ast_corr = False
        if result[1] is True:
            self.ast_corr = True
        
        final_result = result[0]
        
        if not is_parallel:
            if do_norm:
                filter2 = self.filter2
                filter1 = self.filter1
                if self.ast_corr:
                    filter2 = '%s_cor' % self.filter2
                    filter1 = '%s_cor' % self.filter1
                fdict = write_results(result[0], self.agb_mod, self.target,
                                      self.outfile_loc, filter2, filter1,
                                      extra_str=self.extra_str)
                [self.__setattr__(k, v) for k, v in fdict.items()]
            final_result = self.__dict__
        
        return final_result
        """
        return

def load_trilegal_catalog(trilegal_output, filter1, filter2,
                          only_keys='default'):
    '''read the trilegal cat.'''
    if only_keys == 'default':
        only_keys = ['logAge','MH','m_ini','logL','logTe','mM0','Av','mbol',
                     'Mact', 'stage', filter2, filter1, '%s_cor' % filter1,
                     '%s_cor' % filter2, 'CO']
    elif type(only_keys) != list:
        only_keys = None

    sgal = rsp.SimGalaxy(trilegal_output, filter1=filter1, filter2=filter2,
                         only_keys=only_keys)
    return sgal


def gather_results(sgal, norm, target, filter1, filter2,
                   mass_met=True, tpagb_lf=True, narratio_dict=None):
    '''gather results into strings or lists of strings for writing.'''
    result_dict = {}

    if tpagb_lf:
        result_dict['lf_line'] = '# norm mag2 mag1 rgb agb srgb sagb inorm\n' + \
            '\n'.join(['%.4f' % norm,
                       ' '.join(['%g' % m for m in sgal.data[filter2]]),
                       ' '.join(['%g' % m for m in sgal.data[filter1]]),
                       ' '.join(['%i' % m for m in narratio_dict['rgb']]),
                       ' '.join(['%i' % m for m in narratio_dict['agb']]),
                       ' '.join(['%i' % m for m in narratio_dict['srgb']]),
                       ' '.join(['%i' % m for m in narratio_dict['sagb']]),
                       ' '.join(['%i' % m for m in narratio_dict['inorm']])])

    if mass_met is True:
        sgal.all_stages('TPAGB')
        inds = sgal.itpagb
        mag = sgal.mag2

        key = 'mass_met_line'
        mass = sgal.data.m_ini[inds]
        mag = mag[inds]
        mh = sgal.data.MH[inds]
        result_dict[key] = \
            '\n'.join([' '.join(['%g' % t for t in mag]),
                       ' '.join(['%g' % t for t in mass]),
                       ' '.join(['%.3f' % t for t in mh])])

    # N agb/rgb ratio file
    if len(narratio_dict) > 0:
        narratio_fmt = '%(target)s %(nrgb)i %(nagb)i '
        narratio_fmt += '%(ar_ratio).3f %(ar_ratio_err).3f'

        rgb = narratio_dict['rgb']
        agb = narratio_dict['agb']

        nrgb = float(len(rgb))
        nagb = float(len(agb))
        out_dict = {'target': target,
                    'ar_ratio': nagb / nrgb,
                    'ar_ratio_err': rsp.utils.count_uncert_ratio(nagb, nrgb),
                    'nrgb': nrgb,
                    'nagb': nagb}
        result_dict['narratio_line'] = narratio_fmt % out_dict

    return result_dict


def write_results(res_dict, agb_mod, target, outfile_loc, filter2, filter1,
                  extra_str=''):
    '''
    Write results of VSFH output dict to files.

    Paramaters
    ----------
    res_dict : dict
        output of run_once keys with %s_line will be written to a file

    agb_mod, target, filter2, extra_str : strings
        file name formatting stings

    outfile_loc : string
        path to write output file

    Returns
    -------
    fdict : dictionary
        file and path to file
        ex: lf_file: <path_to_lf_file>
    '''
    fmt = '%s_%s_%s_%s_%s%s.dat'
    narratio_header = '# target nrgb nagb ar_ratio ar_ratio_err \n'
    fdict = {}
    for key, line in res_dict.items():
        name = key.replace('_line', '')
        fname = (fmt % (agb_mod, target, filter1, filter2, name, extra_str)).lower()
        fname = os.path.join(outfile_loc, fname)
        with open(fname, 'a') as fh:
            if 'narratio' in key:
                fh.write(narratio_header)
            if type(line) == str:
                line = [line]
            [fh.write('%s \n' % l) for l in line]
        fdict['%s_file' % name] = fname
    return fdict


def run_once(cmd_input_file=None, galaxy_input=None, triout=None, rmfiles=False,
             dry_run=False, do_norm=True, do_norm_kw={},
             cols=None, filter1=None, filter2=None, target=None,
             fake_file=None):

    cmd_input_file = do_norm_kw.get('cmd_input_file', cmd_input_file)
    galaxy_input = do_norm_kw.get('galaxy_input', galaxy_input)
    filter1 = do_norm_kw.get('filter1', filter1)
    filter2 = do_norm_kw.get('filter2', filter2)
    target = do_norm_kw.get('target', target)
    fake_file = do_norm_kw.get('fake_file', fake_file)
    ast_corr = do_norm_kw.get('ast_corr', False)
    diag_plot = do_norm_kw.get('diag_plot', False)
    triout =  do_norm_kw.get('triout', triout)

    rsp.trilegal.utils.run_trilegal(cmd_input_file, galaxy_input, triout,
                                    rmfiles=rmfiles, dry_run=dry_run)
    rsp.trilegal.utils.trilegal2hdf5(triout)
    """
    if ast_corr is True and dry_run is False:
        assert fake_file is not None, 'Need fake file for ast corrections'
        logger.info('adding ast corrections to {}'.format(triout))
        sgal = load_trilegal_catalog(triout, filter1, filter2, only_keys=None)
        rsp.ast_correct_starpop(sgal, overwrite=True, outfile=triout,
                                fake_file=fake_file, diag_plot=False)
        do_norm_kw['sgal'] = sgal

    if do_norm:
        do_norm_kw['triout'] = triout
        sgal, norm, inorm, (srgb, sagb), (rgb, agb) = \
            do_normalization(**do_norm_kw)
        if ast_corr:
            filter1 = '%s_cor' % filter1
            filter2 = '%s_cor' % filter2

        result_dict = gather_results(sgal, norm, target, filter1, filter2,
                                     narratio_dict={'rgb': rgb, 'agb': agb,
                                                    'srgb': srgb, 'sagb': sagb,
                                                    'inorm': inorm})
        result_dict['contam_line'] = contamination_by_phases(sgal, srgb, sagb,
                                                             filter2)
        #if diag_plot:
            #model_cmd_withasts(sgal, rgb=rgb, agb=agb, inorm=inorm,
            #                   **do_norm_kw)
            #import pdb; pdb.set_trace()
        return result_dict, ast_corr
    """

def contamination_by_phases(sgal, srgb, sagb, filter2, diag_plot=False,
                            color_cut=None, target='', line=''):

    """
    contamination by other phases than rgb and agb
    """
    regions = ['MS', 'RGB', 'HEB', 'BHEB', 'RHEB', 'EAGB', 'TPAGB']
    if line == '':
        line += '# %s %s \n' % (' '.join(regions), 'Total')

    sgal.all_stages()
    indss = [sgal.__getattribute__('i%s' % r.lower()) for r in regions]
    try:
        if np.sum(indss) == 0:
            msg = 'No stages in StarPop. Run trilegal with -l flag'
            logger.warning(msg)
            return '{}\n'.format(msg)
    except:
        pass
    if diag_plot is True:
        fig, ax = plt.subplots()

    if color_cut is None:
        inds = np.arange(len(sgal.data[filter2]))
    else:
        inds = color_cut
    mag = sgal.data[filter2][inds]
    ncontam_rgb = [list(set(s) & set(inds) & set(srgb)) for s in indss]
    ncontam_agb = [list(set(s) & set(inds) & set(sagb)) for s in indss]

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
        [ax.hist(mags, bins=bins, alpha=0.5, stacked=True,
                     label=regions)]

    nrgb_cont = np.array([len(n) for n in ncontam_rgb], dtype=int)
    nagb_cont = np.array([len(n) for n in ncontam_agb], dtype=int)

    line += 'rgb %s %i \n' % ( ' '.join(map(str, nrgb_cont)), np.sum(nrgb_cont))
    line += 'agb %s %i \n' % (' '.join(map(str, nagb_cont)), np.sum(nagb_cont))

    line += '# rgb eagb contamination: %i \n' % rheb_eagb_contam
    line += '# frac of total in rgb region: %.3f \n' % frac_rheb_eagb
    line += '# rc contamination: %i \n' % heb_rgb_contam
    line += '# frac of total in rgb region: %.3f \n' % frac_heb_rgb_contam

    logger.info(line)
    if diag_plot is True:
        ax.legend(numpoints=1, loc=0)
        ax.set_title(target)
        plt.savefig('contamination_%s.png' % target, dpi=150)
    return line


def do_normalization(filter1=None, filter2=None, sgal=None, triout=None,
                     offset=None, trgb_exclude=None, trgb=None, ast_corr=False,
                     col_min=None, col_max=None, nrgbs=None, **kwargs):
    '''Do the normalization and save small part of outputs.'''
    if sgal is None:
        assert triout is not None, \
            'need sgal loaded or pass trilegal catalog file name'
        sgal = load_trilegal_catalog(triout, filter1, filter2)

    if ast_corr:
        if not filter1.endswith('cor'):
            filter1 += '_cor'
        if not filter2.endswith('cor'):
            filter2 += '_cor'

    # select rgb and agb regions
    mag_bright = kwargs.get('mag_bright')
    mag_faint = kwargs.get('mag_faint')
    srgb, sagb = rgb_agb_regions(offset, trgb_exclude, trgb, sgal.data[filter2],
                                 col_min=col_min, col_max=col_max,
                                 mag_bright=mag_bright, mag_faint=mag_faint,
                                 mag1=sgal.data[filter1])

    # normalization
    norm, inorm, rgb, agb = normalize_simulation(sgal.data[filter2], nrgbs,
                                                 srgb, sagb)
    return  sgal, norm, inorm, (srgb, sagb), (rgb, agb)

