#!/usr/bin/python
"""
prepare input file

sfh_file: hard coded
TRGB: from table or from Mtrgb
Av, Dmod: From match
RGB region: Color limits by hand, TRGB / completeness or offset
Write optical LF
write ratio table
"""
from __future__ import print_function

import argparse
import difflib
import logging
import os
import sys

import ResolvedStellarPops as rsp
import numpy as np

from ResolvedStellarPops.tpagb_path_config import tpagb_path

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)


if __name__ == '__main__':
    from pop_synth.stellar_pops import limiting_mag, rgb_agb_regions
else:
    from .pop_synth.stellar_pops import limiting_mag, rgb_agb_regions

angst_data = rsp.angst_tables.angst_data


def load_sim_masses(target):
    '''
    adapted from thesis spaz.

    the simulation should have more than 2.5 times the number of stars in
    the CMD as are in the data. Set here are object_mass that should give
    at least that number of stars based on the best fit sfh.
    '''
    if target in ['ngc3741', 'eso540-030', 'ugc-4305-1', 'kkh37', 'ugc-4305-2',
                  'ngc404', 'ngc2976-deep', 'ngc4163', 'ddo78', 'ngc2403-deep']:
        mass = 5e+08
    elif target in ['ddo82', 'ic2574-sgs']:
        mass = 2.5e+09
    elif target in ['ugc-5139']:
        mass = 1.0e+09
    else:
        logger.warning('no info on object mass for {}, assuming 1e8Msun'.format(target))
        mass = 1.0e+08
    return mass


def possible_inputs():
    return {'Av': 0.,
            'binary_frac': 0.35,
            'dmod': 0.,
            'file_origin': None,
            'filter1': None,
            'filter2': None,
            'trgb': None,
            'Mtrgb': None,
            'nrgbs': None,
            'offset': None,
            'photsys': 'wfc3snap',
            'target': None,
            'trgb_exclude': .1,
            'col_max': None,
            'comp_mag1': None,
            'comp_mag2': None,
            'mag_bright': None,
            'col_min': None,
            'mag_faint': None,
            'fake_file': None,
            'sfh_file': None,
            'comp_frac': 0.9,
            'object_sfr_file': None,
            'object_mass': None,
            'outfile_loc': os.getcwd(),
            'file_imf': None,
            'object_cutoffmass': None,
            'cmd_input_file': 'cmd_input_parsecCAF09_V1.2S_M36_S12D2.dat',
            'nsfhs': 1}

def prepare_from_directory(args, search_str, inp_extra):
    """
    Make a partial input file culled from information in several different files
    inp_extra should be used to specify filter1 if more than one filter is in
    a directory.
    """
    args.name = os.path.abspath(args.name)
    assert os.path.isdir(args.name), 'Must supply valid directory name'

    # get matchphot, fake, sfh_file, Hybric MC file
    pars = {'matchphot': \
                rsp.fileio.get_files(args.name,
                                     search_str.format('match'))[0],
            'fake_file': rsp.fileio.get_files(args.name,
                                              search_str.format('fake'))[0],
            'sfh_file': rsp.fileio.get_files(args.name,
                                             search_str.format('sfh'))[0]}

    try:
        pars['hmc_file'] = rsp.fileio.get_files(args.name,
                                                search_str.format('zc'))[0]
        pars['file_origin'] = 'match-hmc'
    except:
        pars['hmc_file'] = pars['sfh_file']
        pars['file_origin'] = 'match-grid'

    # get maglimits from match parameter file
    matchpar, = rsp.fileio.get_files(args.name,
                                     search_str.format('param'))
    # assuming mag is mag2 (i.e, I)
    _, pars['mag_limit_val'] = \
        np.array(open(matchpar).readlines()[5].split()[:-1], dtype=float)

    gal_file = os.path.join(tpagb_path,
                            'SNAP/tables/snap_galaxies.dat')
    gal_table = rsp.fileio.readfile(gal_file,
                                    string_column=[0, -6, -5, -4, -3, -2, -1],
                                    string_length=216)

    target = os.path.split(args.name)[1]
    targets = difflib.get_close_matches(target, gal_table['target'])
    itarg = [i for i, t in enumerate(targets) if t == target]
    if len(itarg) == 1:
        target = targets[itarg[0]]
    
    if len(itarg) > 1:
        assert args.filter is not None, \
            'More than one filter found for {} must choose filter'.format(target)
        
    logger.info('using target: {}'.format(target))

    row = gal_table[np.where(gal_table['target']==target)]
    if args.filter is not None:
        row = row[[i for i, r in enumerate(row)
                   if args.filter.lower() in r['opt1']]]    
        logger.info('user supplied {}, using filter: {}'.format(args.filter,
                                                                row['opt1']))
    # write partial varysfh input file
    newdir = os.path.join(tpagb_path, 'SNAP/varysfh', target)
    rsp.fileio.ensure_dir(newdir)

    pars.update({'outfile_loc': newdir, 'col_min': row['colmin'],
                 'col_max': row['colmax']})

    partial_inpfile = os.path.join(newdir, '{0}{1}.inp'.format(target,
                                                               inp_extra))
    inp = rsp.fileio.InputParameters()
    inp.add_params(pars)
    inp.write_params(partial_inpfile)
    return partial_inpfile

def prepare_for_varysfh(inps, outfile):
    """
    Prepare a default varysfh inputfile
    """
    # get completeness mags for mag limits
    if inps.offset is None:
        logger.info('finding completeness fraction from fake file')
        inps.comp_mag1, inps.comp_mag2 = limiting_mag(inps.fake_file,
                                                      inps.comp_frac)

    # get Av, dmod, mTRGB, photsys
    msfh = rsp.match.utils.MatchSFH(inps.sfh_file)
    inps.Av = msfh.Av
    inps.dmod = msfh.dmod
    angst_target = \
        difflib.get_close_matches(inps.target.upper(),
                                  angst_data.targets)[0].replace('-', '_')
    try:
        target_row = angst_data.__getattribute__(angst_target)
        inps.trgb = target_row['%s,%s' % (inps.filter1, inps.filter2)]['mTRGB']
    except AttributeError:
        logger.error('{} not found in angst tables, \
                     using M=-4 to find mTRGB'.format(angst_target))
        inps.trgb = rsp.astronomy_utils.Mag2mag(-4., inps.filter2, inps.photsys,
                                                dmod=inps.dmod, Av=inps.Av)

    # get observered number of RGB and AGB stars
    mag1, mag2 = np.loadtxt(inps.matchphot, unpack=True)

    if inps.mag_faint is None:
        logger.info('assuming faint mag offset of 2')
        inps.offset = 2.
    else:
        inps.offset = inps.trgb - inps.mag_faint

    rgbs, agbs = rgb_agb_regions(inps.offset, inps.trgb_exclude, inps.trgb,
                                 mag2, mag_faint=inps.mag_faint,
                                 mag_bright=inps.mag_bright,
                                 col_min=inps.col_min, mag1=mag1,
                                 col_max=inps.col_max)

    inps.nrgbs = len(rgbs)
    inps.nagbs = len(agbs)

    # convert match SFH to trilegal AMR
    gal_inp = prepare_galaxy_inputfile(inps)

    inps.__dict__.update(gal_inp)
    inps.write_params(outfile)

def prepare_galaxy_inputfile(inps):
    """Make a galaxy input file for trilegal"""
    # If match was run with setz, this is the logz dispersion.
    # Only useful for clusters, also it is not saved in the match output files
    # Only set in the match parameter file.
    inps.match_zdisp = 0.00
    rsp.match.utils.process_match_sfh(inps.sfh_file,
                                      outfile=inps.object_sfr_file,
                                      zdisp=inps.match_zdisp)

    gal_dict =\
        {'mag_limit_val': limiting_mag(inps.fake_file, 0.1)[1],
         'object_av': inps.Av,
         'object_dist': 10 ** (inps.dmod/5. + 1.),
         'photsys': inps.photsys,
         'object_mass': inps.object_mass or load_sim_masses(inps.target),
         'object_sfr_file': inps.object_sfr_file,
         'file_imf': inps.file_imf or 'tab_imf/imf_kroupa_match.dat',
         'binary_frac': inps.binary_frac or 0.35,
         'object_cutoffmass': inps.object_cutoffmass or 0.8}

    # filter1 is used here to find the mag depth for trilegal input.
    gal_dict['filter1'] = inps.filter2
    trigal_dict = rsp.trilegal.utils.galaxy_input_dict(**gal_dict)
    del gal_dict['filter1']

    gal_inp = rsp.fileio.InputParameters(default_dict=trigal_dict)
    gal_inp.write_params(inps.galaxy_input,
                         rsp.trilegal.utils.galaxy_input_fmt())
    
    return gal_dict

def prepare_outfiles(inps, inp_extra):
    """
    set up default output file syntax
    """
    # where things are going
    if inps.outfile_loc is None:
        inps.outfile_loc = os.path.split(inps.matchphot)[0]
        assert os.path.isdir(inps.outfile_loc), 'bad output directory'
    else:
        rsp.fileio.ensure_dir(inps.outfile_loc)

    inps.target, inps.filter1, inps.filter2 = rsp.parse_pipeline(inps.matchphot)
    # trilegal sfr filename
    inps.object_sfr_file = os.path.join(inps.outfile_loc,
                                        '{0}{1}.trisfr'.format(inps.target,
                                                               inp_extra))
    # trilegal galaxy input filename
    inps.galaxy_input = os.path.join(inps.outfile_loc,
                                     '{0}{1}.galinp'.format(inps.target,
                                                            inp_extra))
    return inps

def main(argv):
    """
    need matchphot in the format of target filter1 filter2
    need match sfh file for Av and dmod
    col_min, col_max, mag_faint, mag_bright
    fakefile
    photsys
    'object_mass': inps.object_mass or 1e7,
    'object_sfr_file': inps.object_sfr_file,
    'file_imf': inps.file_imf or 'tab_imf/imf_salpeter.dat',
    'binary_frac': inps.binary_frac or 0.,
    'object_cutoffmass': inps.object_cutoffmass or 0.8}
    """

    parser = argparse.ArgumentParser(description="Create input file \
                                                  for VarySFH")

    parser.add_argument('-d', '--directory', action='store_true',
                        help='create partial input file from \
                              specified directory')

    parser.add_argument('-v', '--pdb', action='store_true',
                        help='verbose mode')

    parser.add_argument('-f', '--filter', type=str, default=None,
                        help='V filter (if more than one in directory)')

    parser.add_argument('-n', '--nsfhs', type=str, default=1,
                        help='Number of sampled SFHs to run')
    
    parser.add_argument('name', type=str,
                        help='partial input file or if using -d, \
                             directory name')

    args = parser.parse_args(argv)

    # set up logging
    handler = logging.FileHandler('{}_prepare_data.log'.format(args.name))
    if args.pdb:
        handler.setLevel(logging.DEBUG)
    else:
        handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    if args.directory:
        if args.filter is not None:
            fsearch = '*{}'.format(args.filter)
            inp_extra = '_{}'.format(args.filter)
        else:
            fsearch = ''
            inp_extra = ''
        search_str = fsearch + '*{}'
        partial_inpfile = prepare_from_directory(args, search_str, inp_extra)
    else:
        partial_inpfile = args.name

    inputs = dict(possible_inputs().items() + {'nsfhs': args.nsfhs}.items())
    inps = rsp.fileio.InputParameters(default_dict=inputs)
    inps.add_params(rsp.fileio.load_input(partial_inpfile), loud=args.pdb)

    inps = prepare_outfiles(inps, inp_extra)

    # varysfh input file name
    outfile = partial_inpfile.replace('.inp', '.vsfhinp')
            
    prepare_for_varysfh(inps, outfile)
    return


if __name__ == "__main__":
    main(sys.argv[1:])
