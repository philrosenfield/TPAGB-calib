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
from ..sfhs.star_formation_histories import StarFormationHistories
from ..TPAGBparams import snap_src


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

def load_trilegal_catalog(trilegal_catalog):
    '''read a table into a SimGalaxy object'''
    return rsp.SimGalaxy(trilegal_catalog)


def make_ast_corrections(trilegal_catalogs, target, outfile='default',
                         diag_plot=False, overwrite=True):
    """
    apply ast corrections from fake files found in matchfake_loc/*[target]*
    see rsp.ast_correct_starpop
    """
    # where the matchfake files live
    matchfake_loc = os.path.join(snap_src, 'data', 'galaxies')
    
    # search string for fake files
    search_str = '*{}*.matchfake'.format(target.upper())
    
    fakes = rsp.fileio.get_files(matchfake_loc, search_str)
    logger.info('fake files found: {}'.format(fakes))
    asts = [rsp.ASTs(f) for f in fakes]
    
    for trilegal_catalog in trilegal_catalogs:
        logger.info('working on {}'.format(trilegal_catalog))
        sgal = load_trilegal_catalog(trilegal_catalog)
        # "overwrite" (append columns) to the existing catalog by default
        if outfile == 'default':
            outfile = trilegal_catalog

        # do the ast corrections
        [rsp.ast_correct_starpop(sgal, asts_obj=ast, overwrite=overwrite,
                                 outfile=outfile, diag_plot=diag_plot)
         for ast in asts]
    return


def do_normalization(filter1=None, filter2=None, sgal=None, tricat=None,
                     offset=None, trgb_exclude=None, trgb=None, nrgbs=None,
                     col_min=None, col_max=None, mag_bright=None,
                     mag_faint=None):
    '''Do the normalization and save small part of outputs.'''
    if sgal is None:
        sgal = load_trilegal_catalog(tricat)

    # select rgb and agb regions
    sgal_rgb, sgal_agb = rgb_agb_regions(offset, trgb_exclude, trgb,
                                         sgal.data[filter2], col_min=col_min,
                                         col_max=col_max, mag_bright=mag_bright,
                                         mag_faint=mag_faint,
                                         mag1=sgal.data[filter1])

    # normalization
    norm, idx_norm, sim_rgb, sim_agb = normalize_simulation(sgal.data[filter2],
                                                            nrgbs, sgal_rgb,
                                                            sgal_agb)

    return sgal, norm, idx_norm, (sgal_rgb, sgal_agb), (sim_rgb, sim_agb)



# main
#-d directory or trilegal catalog name
# just do ast_correction and leave
# output lf or cmd optical or ir
# settings for normalization
# read in observations table

def main(argv):
    # indict:
    # ast_corr bool
    # fake_file if ast_corr true
    # filter1
    # filter2
    # target
    # agb_mod
    # outfile loc
    # extra_str
    # nrgbs
    # offset,
    # trgb_exclude,
    # trgb,
    # col_min,
    # col_max,
    # mag_bright,
    # mag_faint,
    # need ir stuff too! -- fill out table...
    parser = argparse.ArgumentParser(description="Cull useful info from \
                                                  trilegal catalog")

    parser.add_argument('-d', '--directory', action='store_true',
                        help='opperate on all files in a directory')

    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose mode')

    parser.add_argument('-a', '--ast', action='store_true',
                        help='make ast corrections to file(s)')

    parser.add_argument('-t', '--target', type=str, help='target name')

    parser.add_argument('name', type=str, nargs='*',
                        help='trilegal catalog or directory if -d flag')

    args = parser.parse_args(argv)

    # set up logging
    handler = logging.FileHandler('analyze.log')
    if args.verbose:
        handler.setLevel(logging.DEBUG)
    else:
        handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    import pdb; pdb.set_trace()
    if args.directory:
        tricats = rsp.fileio.get_files(args.name[0], '*dat')
    else:
        tricats = args.name
    
    if not args.target:   
        if args.directory:
            target = os.path.split(args.name[0])[1]
        else:
            target = tricat.split('_')[1]
        
    if args.verbose:
        logger.info('working on target: {}'.format(target))
        
    make_ast_corrections(tricats, target)

    """
    if indict['ast_corr']:
        if not hasattr(sgal.data, '{}_cor'.format(filter2)):
            assert indict['fake_file'] is not None, \
                'No ast correction in {}, need a fake file to make them on the fly'.format(sgal.name)
            ast_correction(sgal, indict['fake_file'], outfile='default',
                           diag_plot=False, overwrite=True)
            indict['filter1'] += '_cor'
            indict['filter2'] += '_cor'

    sgal, norm, idx_norm, (sgal_rgb, sgal_agb), (sim_rgb, sim_agb) = \
        do_normalization(**indict)
        
    result_dict = gather_results(sgal, norm, indict['target'],
                                 indict['filter1'],
                                 indict['filter2'],
                                 narratio_dict={'sim_rgb': rgb,
                                                'sim_agb': agb,
                                                'sgal_rgb': srgb,
                                                'sgal_agb': sagb,
                                                'idx_norm': inorm})
    
    result_dict['contam_line'] = contamination_by_phases(sgal, sgal_rgb,
                                                         sgal_agb, filter2)

    file_dict = write_results(result_dict, agb_mod, indict['target'],
                              indict['outfile_loc'], indict['filter1'],
                              indict['filter2'])
    
    # add file_dict to plotinp
    """
    return

if __name__ == "__main__":
    main(sys.argv[1:])

def gather_results(sgal, norm, target, filter1, filter2,
                   mass_met=True, tpagb_lf=True, narratio_dict=None):
    '''gather results into strings or lists of strings for writing.'''
    result_dict = {}

    if tpagb_lf:
        result_dict['lf_line'] = '# norm mag2 mag1 sim_rgb sim_agb sgal_rgb sgal_agb idx_norm\n' + \
            '\n'.join(['%.4f' % norm,
                       ' '.join(['%g' % m for m in sgal.data[filter2]]),
                       ' '.join(['%g' % m for m in sgal.data[filter1]]),
                       ' '.join(['%i' % m for m in narratio_dict['rgb']]),
                       ' '.join(['%i' % m for m in narratio_dict['agb']]),
                       ' '.join(['%i' % m for m in narratio_dict['srgb']]),
                       ' '.join(['%i' % m for m in narratio_dict['sagb']]),
                       ' '.join(['%i' % m for m in narratio_dict['inorm']])])

    #if mass_met is True:
    #    sgal.all_stages('TPAGB')
    #    inds = sgal.itpagb
    #    mag = sgal.mag2

    #    key = 'mass_met_line'
    #    mass = sgal.data.m_ini[inds]
    #    mag = mag[inds]
    #    mh = sgal.data.MH[inds]
    #    result_dict[key] = \
    #        '\n'.join([' '.join(['%g' % t for t in mag]),
    #                   ' '.join(['%g' % t for t in mass]),
    #                   ' '.join(['%.3f' % t for t in mh])])

    # N agb/rgb ratio file
    if len(narratio_dict) > 0:
        narratio_fmt = '%(target)s %(filter2)s %(nrgb)i %(nagb)i '
        narratio_fmt += '%(ar_ratio).3f %(ar_ratio_err).3f'

        rgb = narratio_dict['rgb']
        agb = narratio_dict['agb']

        nrgb = float(len(rgb))
        nagb = float(len(agb))
        out_dict = {'target': target,
                    'filter2': filter2,
                    'ar_ratio': nagb / nrgb,
                    'ar_ratio_err': rsp.utils.count_uncert_ratio(nagb, nrgb),
                    'nrgb': nrgb,
                    'nagb': nagb}
        result_dict['narratio_line'] = narratio_fmt % out_dict

    return result_dict

def write_results(res_dict, agb_mod, target, outfile_loc, filter1, filter2,
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


### Snippets below ###

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

