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
import os
import sys

import ResolvedStellarPops as rsp
import numpy as np

from ResolvedStellarPops.tpagb_path_config import tpagb_path

if __name__ == '__main__':
    from pop_synth.stellar_pops import limiting_mag, rgb_agb_regions
else:
    from .pop_synth.stellar_pops import limiting_mag, rgb_agb_regions


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
            'photsys': None,
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
            'object_cutoffmass': None}


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
    angst_data = rsp.angst_tables.angst_data

    parser = argparse.ArgumentParser(description="Create input file for VarySFH")

    parser.add_argument('-d', '--directory', action='store_true',
                        help='create partial input file from specified directory')

    parser.add_argument('-v', '--pdb', action='store_true',
                        help='debugging mode')

    parser.add_argument('-f', '--filter', type=str, default=None
                        help='V filter (if more than one in directory')

    parser.add_argument('name', type=str,
                        help='partial input file or if using -d, directory name')

    args = parser.parse_args(argv)

    if args.pdb:
        import pdb
        pdb.set_trace()

    if args.f is not None:
        fsearch = '*{}'.format(args.f)
    else:
        fsearch = ''

    search_str = fsearch + '*{}'
    if args.directory:
        assert os.path.isdir(args.name), 'Must supply valid directory name'
        pars = {'matchphot': rsp.fileio.get_files(args.name, search_str.format('match'))[0],
                'fake_file': rsp.fileio.get_files(args.name, search_str.format('fake'))[0],
                'sfh_file': rsp.fileio.get_files(args.name, search_str.format('sfh'))[0]}

        try:
            pars['hmc_file'] = rsp.fileio.get_files(args.name, search_str.format('zc'))[0]
            pars['file_origin'] = 'match-hmc'
        except:
            pars['hmc_file'] = pars['sfh_file']
            pars['file_origin'] = 'match-grid'

        matchpar, = rsp.fileio.get_files(args.name, search_str.format('param'))
        # assuming mag is mag2 (i.e, I)
        pars['mag_bright'], pars['mag_faint'] = \
            np.array(open(matchpar).readlines()[5].split()[:-1], dtype=float)

        # find matchphot, fake, sfh_file, col_min, col_max, mag_faint,
        # mag_bright write to file
        target = os.path.split(args.name)[1]
        newdir = os.path.join(tpagb_path, 'SNAP/varysfh', target)
        gal_file = os.path.join(tpagb_path,
                                'SNAP/tables/paperII_varsfh_table.dat')
        gal_table = rsp.fileio.readfile(gal_file, string_column=[0, -2, -1],
                                        string_length=216)
        target = difflib.get_close_matches(target, gal_table['target'])[0]
        print('using target: {}'.format(target))
        row = gal_table[np.where(gal_table['target']==target)]
        pars.update({'outfile_loc': newdir, 'col_min': row['colmin'],
                     'col_max': row['colmax']})
        rsp.fileio.ensure_dir(newdir)
        filename = os.path.join(newdir, '%s.inp' % target)
        inp = rsp.fileio.InputParameters()
        inp.add_params(pars)
        inp.write_params(filename)
    else:
        filename = args.name

    inps = rsp.fileio.InputParameters(default_dict=possible_inputs())
    inps.add_params(rsp.fileio.load_input(filename))
    inps.target, inps.filter1, inps.filter2 = rsp.parse_pipeline(inps.matchphot)

    if inps.offset is None:
        print('finding completeness fraction from fake file')
        inps.comp_mag1, inps.comp_mag2 = limiting_mag(inps.fake_file,
                                                      inps.comp_frac)

    msfh = rsp.match.utils.MatchSFH(inps.sfh_file)
    inps.Av = msfh.Av
    inps.dmod = msfh.dmod
    angst_target = \
        difflib.get_close_matches(inps.target.upper(),
                                  angst_data.targets)[0].replace('-', '_')
    try:
        target_row = angst_data.__getattribute__(angst_target)
        inps.trgb = target_row['%s,%s' % (inps.filter1, inps.filter2)]['mTRGB']
        if inps.photsys is None:
            inps.photsys = angst_data.get_item(angst_target, 'camera').lower()
            if inps.photsys == 'acs':
                inps.photsys = 'acs_wfc'
    except AttributeError:
        print('{} not found in angst tables, \
              using M=-4 to find mTRGB'.format(angst_target))
        inps.trgb = rsp.astronomy_utils.Mag2mag(-4., inps.filter2, inps.photsys,
                                                dmod=inps.dmod, Av=inps.Av)

    mag1, mag2 = np.loadtxt(inps.matchphot, unpack=True)

    if inps.mag_faint is None:
        print('assuming faint mag offset of 2')
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
    if inps.outfile_loc is None:
        inps.outfile_loc = os.path.split(inps.matchphot)[0]
    else:
        rsp.fileio.ensure_dir(inps.outfile_loc)

    inps.object_sfr_file = os.path.join(inps.outfile_loc,
                                        '%s.sfr' % inps.target.lower())

    # If match was run with setz, this is the logz dispersion.
    # Only useful for clusters, also it is not saved in the match output files
    # only set in the match parameter file.
    inps.match_zdisp = 0.00
    rsp.match.utils.process_match_sfh(inps.sfh_file,
                                      outfile=inps.object_sfr_file,
                                      zdisp=inps.match_zdisp)
    inps.galaxy_input = os.path.join(inps.outfile_loc,
                                     '%s.inp' % inps.target.lower())

    gal_inp =\
        {'mag_limit_val': limiting_mag(inps.fake_file, 0.1)[1],
         'object_av': inps.Av,
         'object_dist': 10 ** (inps.dmod/5. + 1.),
         'photsys': inps.photsys,
         'object_mass': inps.object_mass or 1e7,
         'object_sfr_file': inps.object_sfr_file,
         'file_imf': inps.file_imf or 'tab_imf/imf_salpeter.dat',
         'binary_frac': inps.binary_frac or 0.,
         'object_cutoffmass': inps.object_cutoffmass or 0.8}

    inps.__dict__.update(gal_inp)
    gal_inp['filter1'] = inps.filter2
    gal_inp = rsp.trilegal.utils.galaxy_input_dict(**gal_inp)

    gal_inp = rsp.fileio.InputParameters(default_dict=gal_inp)
    gal_inp.write_params(inps.galaxy_input,
                         rsp.trilegal.utils.galaxy_input_fmt())
    inps.write_params(filename.replace('.inp', 'prepped.inp'))
    return inps

if __name__ == "__main__":
    main(sys.argv[1:])
