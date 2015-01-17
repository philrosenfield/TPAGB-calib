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
if __name__ == '__main__':
    from pop_synth.stellar_pops import limiting_mag, rgb_agb_regions
else:
    from .pop_synth.stellar_pops import limiting_mag, rgb_agb_regions
import ResolvedStellarPops as rsp
import numpy as np
import sys
import os
import argparse


def possible_inputs():
    return {'Av': 0.,
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
    
    parser.add_argument('-d', '--directory', type=str, action='store_true',
                        help='act on a MATCH directory instead of a partial input file')
    
    parser.add_argument('-v', '--pdb', action='store_true',
                        help='debugging mode')

    parser.add_argument('name', type=str,
                        help='filename or if using -d, directory name')
    
    args = parser.parse_args(argv)
     
    if args.v:
        import pdb; pdb.set_trace()
    
    if args.d:
        pars = {'matchphot': rsp.fileio.get_files(args.name, '*match')[0],
                'fake_file': rsp.fileio.get_files(args.name, '*fake')[0],
                'sfh_file': rsp.fileio.get_files(args.name, '*sfh')[0]}
        
        try:
            pars['hmc_file'] = rsp.fileio.get_files(args.name, '*zc')[0]
        except:
            pass
        
        matchpar, = rsp.fileio.get_files(args.name, '*param')
        # assuming mag is mag2 (i.e, I)
        pars['mag_bright'], pars['mag_faint'] = \
            map(float, open(matchpar).readlines()[5].split()[:-1])

        # find matchphot, fake, sfh_file, col_min, col_max, mag_faint, mag_bright
        # write to file
        target = os.path.split(args.name)[1]
        vsfh_loc = '/home/rosenfield/research/TP-AGBcalib/SNAP/varysfh/'
        newdir = os.path.join(vsfh_loc, target)
        rsp.fileio.ensure_dir(newdir)
        filename = os.path.join(newdir, '%s.inp' % target)
        inp = rsp.fileio.InputParameters()
        inp.write_params(pars)
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
    try:
        target_row = angst_data.__getattribute__(inps.target.upper())
        inps.trgb = target_row['%s,%s' % (inps.filter1, inps.filter2)]['mTRGB']
        if inps.photsys is None:
            inps.photsys = angst_data.get_item(inps.target.upper(), 'camera')
    except AttributeError:
        print('{} not found in angst tables, using M=-4 to find mTRGB'.format(inps.target.upper()))
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
    #inps.match_zdisp = 0.05
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
    inp_obj.write_params(filename.replace('.inp', 'prepped.inp'))
    return inps

if __name__ == "__main__":
    main(sys.argv[1:])
