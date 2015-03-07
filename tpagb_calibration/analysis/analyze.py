"""
Utilities to analyze trilegal output catalogs of TPAGB models
"""
import argparse
import logging
import numpy as np
import os
import sys
import time

#import matplotlib as mpl
#mpl.use('Agg')
import matplotlib.pylab as plt
import ResolvedStellarPops as rsp

from IPython import parallel
from ..pop_synth.stellar_pops import normalize_simulation, rgb_agb_regions
from ..plotting.plotting import model_cmd_withasts
from ..sfhs.star_formation_histories import StarFormationHistories
from ..TPAGBparams import snap_src

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

optfilter2 = 'F814W'
nirfilter2 = 'F160W'
nirfilter1 = 'F110W'

data_loc = os.path.join(snap_src, 'data', 'galaxies')

def load_trilegal_catalog(trilegal_catalog):
    '''read a table into a SimGalaxy object'''
    return rsp.SimGalaxy(trilegal_catalog)


def cutheb(sgal):
    """
    Mark HeB stars as unrecovered from observations -- assumes filters start
    with "F" as in F110W or F814W_cor.
    
    i.e.:
    assign 99. to sgal.data[f] where sgal.data['stage'] is between 4 and 6.
    where f is any array that begins with "F"
    
    4-5 is based on stages ordered (indexed 0) =
    PMS, MS, SUBGIANT, RGB, HEB, RHEB, BHEB, EAGB, TPAGB, POSTAGB, WD

    To do:
    Add mbol or logl limit on where to apply this correction
    """
    # HeB, RHeB, BHeB
    sgal.iheb, = np.nonzero((sgal.data['stage'] >= 4) & \
                            (sgal.data['stage'] <= 6))
    
    filters = [f for f in sgal.data.dtype.names if f.startswith('F')]
    
    # flag heb star magnitudes as "unrecovered"
    for f in filters:
        sgal.data[f][sgal.iheb] = 99.
    
    return sgal


def check_astcor(filters):
    """add _cor to filter names if it isn't already there"""
    if type(filters) is str:
        filters = [filters]

    for i, f in enumerate(filters):
        if not f.endswith('cor'):
            filters[i] = f + '_cor'
    return filters


def do_normalization(opt=True, ast_cor=True, optfilter1=None, sgal=None,
                     tricat=None, nrgbs=None, cut_heb=False, regions_kw={}):
    '''Do the normalization: call rgb_agb_regions and normalize_simulations.'''

    if opt:
        filter1 = optfilter1
        filter2 = optfilter2
    else:
        filter1 = nirfilter1
        filter2 = nirfilter2

    if ast_cor:
        filter1, filter2 = check_astcor([filter1, filter2])

    if sgal is None:
        sgal = load_trilegal_catalog(tricat)

    if cut_heb:
        sgal = cutheb(sgal)

    # select rgb and agb regions
    sgal_rgb, sgal_agb = rgb_agb_regions(sgal.data[filter2],
                                         mag1=sgal.data[filter1],
                                         **regions_kw)

    # normalization
    norm, idx_norm, sim_rgb, sim_agb = normalize_simulation(sgal.data[filter2],
                                                            nrgbs, sgal_rgb,
                                                            sgal_agb)

    if opt:
        norm_dict = {'optnorm': norm,
                     'optsim_rgb': sim_rgb,
                     'optsim_agb': sim_agb,
                     'optsgal_rgb': sgal_rgb,
                     'optsgal_agb': sgal_agb,
                     'optidx_norm': idx_norm}
    else:
        norm_dict = {'nirnorm': norm,
                     'nirsim_rgb': sim_rgb,
                     'nirsim_agb': sim_agb,
                     'nirsgal_rgb': sgal_rgb,
                     'nirsgal_agb': sgal_agb,
                     'niridx_norm': idx_norm}
        
    return sgal, norm_dict


def makelf(trilegal_catalogs, target, heb=True, norm=True, ast_cor=True,
           completeness=True, data=True, norm_kw={}, lf_line=''):
    pass


def tpagb_lf(sgal, narratio_dict, optfilt1, optfilt2, nirfilt1, nirfilt2,
             lf_line=''):
    """format a narratio_dict for a line in the LF output file"""
    optrgb = narratio_dict['optsim_rgb']
    optagb = narratio_dict['optsim_agb']
    nirrgb = narratio_dict['nirsim_rgb']
    niragb = narratio_dict['optsim_agb']

    header = '# {} {} {} {} '.format(optfilt1, optfilt2, nirfilt1,
                                     nirfilt2)
    header += 'optsim_rgb optsim_agb optsgal_rgb optsgal_agb '
    header += 'optidx_norm optnorm nirsim_rgb nirsim_agb nirsgal_rgb '
    header += 'nirsgal_agb niridx_norm nirnorm\n'
    
    if len(lf_line) == 0:
        lf_line = header
    lf_line += '\n'.join([' '.join(['%g' % m for m in sgal.data[optfilt1]]),
                          ' '.join(['%g' % m for m in sgal.data[optfilt2]]),
                          ' '.join(['%g' % m for m in sgal.data[nirfilt1]]),
                          ' '.join(['%g' % m for m in sgal.data[nirfilt2]]),
                          ' '.join(['%i' % m for m in optrgb]),
                          ' '.join(['%i' % m for m in optagb]),
                          ' '.join(['%i' % m for m in narratio_dict['optsgal_rgb']]),
                          ' '.join(['%i' % m for m in narratio_dict['optsgal_agb']]),
                          ' '.join(['%i' % m for m in narratio_dict['optidx_norm']]),
                          '%.4f' % narratio_dict['optnorm'],
                          ' '.join(['%i' % m for m in nirrgb]),
                          ' '.join(['%i' % m for m in niragb]),
                          ' '.join(['%i' % m for m in narratio_dict['nirsgal_rgb']]),
                          ' '.join(['%i' % m for m in narratio_dict['nirsgal_rgb']]),
                          ' '.join(['%i' % m for m in narratio_dict['niridx_norm']]),
                          '%.4f' % narratio_dict['nirnorm']])
    return lf_line


def narratio(target, optnrgb, optnagb, nirnrgb, nirnagb, optfilt2, nirfilt2,
             narratio_line=''):
    """format numbers of stars for the narratio table"""
    # N agb/rgb ratio file
    narratio_fmt = '%(target)s %(optfilter2)s %(optnrgb)i %(optnagb)i '
    narratio_fmt += '%(optar_ratio).3f %(optar_ratio_err).3f '
    narratio_fmt += '%(nirfilter2)s %(nirnrgb)i %(nirnagb)i '
    narratio_fmt += '%(nirar_ratio).3f %(nirar_ratio_err).3f\n'

    out_dict = {'target': target,
                'optfilter2': optfilt2,
                'optar_ratio': optnagb / optnrgb,
                'optar_ratio_err': rsp.utils.count_uncert_ratio(optnagb, optnrgb),
                'optnrgb': optnrgb,
                'optnagb': optnagb,
                'nirfilter2': nirfilt2,
                'nirar_ratio': nirnagb / nirnrgb,
                'nirar_ratio_err': rsp.utils.count_uncert_ratio(nirnagb, nirnrgb),
                'nirnrgb': nirnrgb,
                'nirnagb': nirnagb}
    narratio_line += narratio_fmt % out_dict
    return narratio_line


def gather_results(sgal, target, optfilter1, ast_cor=True, narratio_dict=None,
                   lf_line='', narratio_line=''):
    '''gather results into strings: call tpagb_lf and narratio'''
    if ast_cor:
        optfilt1, optfilt2, nirfilt1, nirfilt2 = check_astcor([optfilter1,
                                                               optfilter2,
                                                               nirfilter1,
                                                               nirfilter2])

    lf_line = tpagb_lf(sgal, narratio_dict, optfilt1, optfilt2, nirfilt1, nirfilt2,
                       lf_line=lf_line)

    optrgb = narratio_dict['optsim_rgb']
    optagb = narratio_dict['optsim_agb']
    nirrgb = narratio_dict['nirsim_rgb']
    niragb = narratio_dict['optsim_agb']
    
    optnrgb = float(len(optrgb))
    optnagb = float(len(optagb))
    nirnrgb = float(len(nirrgb))
    nirnagb = float(len(niragb))
    
    narratio_line = narratio(target, optnrgb, optnagb, nirnrgb, nirnagb,
                             optfilt2, nirfilt2, narratio_line=narratio_line)

    return lf_line, narratio_line


def write_results(res_dict, target, outfile_loc, optfilter1, extra_str=''):
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
    narratio_header = '# target optnrgb optnagb optar_ratio optar_ratio_err '
    narratio_header += 'nirnrgb nirnagb nirar_ratio nirar_ratio_err \n'

    fdict = {}
    for key, line in res_dict.items():
        name = key.replace('_line', '')
        fname = ('_'.join(['%s' % s for s in (target, optfilter1,
                                              optfilter2, nirfilter1,
                                              nirfilter2, name)])).lower()

        fname = os.path.join(outfile_loc, '%s%s.dat' % (fname, extra_str))
        with open(fname, 'a') as fh:
            if 'narratio' in key:
                fh.write(narratio_header)
            if type(line) == str:
                line = [line]
            [fh.write('%s \n' % l) for l in line]
        fdict['%s_file' % name] = fname
    return fdict


def get_trgb(target, optfilter1=None):
    import difflib
    angst_data = rsp.angst_tables.angst_data

    angst_target = difflib.get_close_matches(target.upper(),
                                             angst_data.targets)[0].replace('-', '_')
    
    target_row = angst_data.__getattribute__(angst_target)
    opt_trgb = target_row['%s,%s' % (optfilter1, optfilter2)]['mTRGB']

    nir_trgb = angst_data.get_snap_trgb_av_dmod(target.upper())[0]
    return opt_trgb, nir_trgb


def laad_obs(target, optfilter1=''):
    """load in NIR and OPT galaxy as StarPop objects"""
    from astropy.io import fits
    nirgalname, = rsp.fileio.get_files(data_loc, '*{}*fits'.format(target.upper()))
    optgalname, = rsp.fileio.get_files(data_loc, ('*{}*{}*fits'.format(target, optfilter1).lower()))
    nirgal = rsp.StarPop()
    nirgal.data = fits.getdata(nirgalname)

    optgal = rsp.StarPop()
    optgal.data = fits.getdata(optgalname)
    return optgal, nirgal


def main(argv):
    parser = argparse.ArgumentParser(description="Cull useful info from \
                                                  trilegal catalog")

    parser.add_argument('-a', '--ast_cor', action='store_true',
                        help='use ast corrected mags')
    
    parser.add_argument('-c', '--colorlimits', type=str, default=None,
                        help='comma separated color min, color max, opt then nir')

    parser.add_argument('-d', '--directory', action='store_true',
                        help='opperate on *hdf5 files in a directory')

    parser.add_argument('-e', '--trgbexclude', type=str, default='0.1,0.2',
                        help='comma separated regions around trgb to exclude')

    parser.add_argument('-f', '--optfilter1', type=str,
                        help='optical V filter')

    parser.add_argument('-z', '--cut_heb', action='store_true',
                        help='cut HeB files from analysis')

    parser.add_argument('-m', '--maglimits', type=str, default=None,
                        help='comma separated mag faint, mag bright, opt then nir')

    parser.add_argument('-o', '--trgboffsets', type=str, default='2.0,1.5',
                        help='comma separated trgb offsets')

    parser.add_argument('-t', '--target', type=str, help='target name')
    
    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose mode')
    
    parser.add_argument('name', type=str, nargs='*',
                        help='trilegal catalog(s) or directory if -d flag')

    args = parser.parse_args(argv)


    if not args.target:   
        if args.directory:
            target = os.path.split(args.name[0])[1]
        else:
            target = args.name[0].split('_')[1]
    else:
        target = args.target

    if args.directory:
        tricats = rsp.fileio.get_files(args.name[0], '*_???.hdf5')
        outfile_loc = args.name[0]
    else:
        tricats = args.name
        outfile_loc = os.path.split(args.name[0])[0]
        
    # set up logging
    logfile = os.path.join(outfile_loc, '{}_analyze.log'.format(target))
    handler = logging.FileHandler(logfile)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    if args.verbose:
        logger.setLevel(logging.DEBUG)
    else:
        logger.setLevel(logging.INFO)    
    
    logger.info('command: {}'.format(' '.join(argv)))
    logger.info('logfile: {}'.format(logfile))
    logger.debug('working on target: {}'.format(target))

    # need the following in opt and nir
    if args.trgboffsets is not None:
        opt_offset, nir_offset = map(float, args.trgboffsets.split(','))
    else:
        opt_offset, nir_offset = None, None
    
    if args.colorlimits is not None:
        opt_colmin, opt_colmax, nir_colmin, nir_colmax = \
            map(float, args.colorlimits.split(','))
    else:
        opt_colmin, opt_colmax, nir_colmin, nir_colmax = None, None, None, None

    if args.maglimits is not None:
        opt_magfaint, opt_magbright, nir_magfaint, nir_magbright = \
            map(float, args.colorlimits.split(','))
    else:
        opt_magfaint, opt_magbright, nir_magfaint, nir_magbright = None, None, None, None

    opt_trgbexclude, nir_trgbexclude = map(float, args.trgbexclude.split(','))

    # get the number of rgb and agb stars from the observations
    opt_trgb, nir_trgb = get_trgb(target, optfilter1=args.optfilter1)
    optgal, nirgal = laad_obs(target, optfilter1=args.optfilter1)

    optregions_kw = {'offset': opt_offset,
                     'trgb_exclude': opt_trgbexclude,
                     'trgb': opt_trgb,
                     'col_min': opt_colmin,
                     'col_max': opt_colmax,
                     'mag_bright': opt_magbright,
                     'mag_faint': opt_magfaint}

    nirregions_kw = {'offset': nir_offset,
                     'trgb_exclude': nir_trgbexclude,
                     'trgb': nir_trgb,
                     'col_min': nir_colmin,
                     'col_max': nir_colmax,
                     'mag_bright': nir_magbright,
                     'mag_faint': nir_magfaint}

    optgal_rgb, optgal_agb = rgb_agb_regions(optgal.data['MAG2_ACS'],
                                             mag1=optgal.data['MAG1_ACS'],
                                             **optregions_kw)

    nirgal_rgb, nirgal_agb = rgb_agb_regions(nirgal.data['MAG2_IR'],
                                             mag1=nirgal.data['MAG1_IR'],
                                             **nirregions_kw)
    
    obs_optnrgbs = float(len(optgal_rgb))
    obs_nirnrgbs = float(len(nirgal_rgb))
    obs_optnagbs = float(len(optgal_agb))
    obs_nirnagbs = float(len(nirgal_agb))

    narratio_line = narratio('data', obs_optnrgbs, obs_optnagbs, obs_nirnrgbs,
                             obs_nirnagbs, optfilter2, nirfilter2)

    # normalize each trilegal catalog    
    norm_kws = {'cut_heb': args.cut_heb, 'ast_cor': args.ast_cor}
    lf_line = ''
    for tricat in tricats:
        logger.debug('normalizing: {}'.format(tricat))
        sgal, optnorm_dict = do_normalization(opt=True, tricat=tricat,
                                              optfilter1=args.optfilter1,
                                              nrgbs=obs_optnrgbs,
                                              regions_kw=optregions_kw,
                                              **norm_kws)

        sgal, nirnorm_dict = do_normalization(opt=False, sgal=sgal,
                                              nrgbs=obs_nirnrgbs,
                                              regions_kw=nirregions_kw,
                                              **norm_kws)

        narratio_dict = dict(optnorm_dict.items() + nirnorm_dict.items())
        
        lf_line, narratio_line = gather_results(sgal, target, args.optfilter1,
                                                narratio_dict=narratio_dict,
                                                lf_line=lf_line,
                                                narratio_line=narratio_line)
    
    result_dict = {'lf_line': lf_line, 'narratio_line': narratio_line}
    #result_dict['contam_line'] = contamination_by_phases(sgal, sgal_rgb,
    #                                                     sgal_agb, filter2)
    
    # write the output files
    file_dict = write_results(result_dict, target, outfile_loc,
                              args.optfilter1)
    print file_dict
    # plot?
    return


if __name__ == "__main__":
    main(sys.argv[1:])


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

