"""
prepare input file

sfh_file: hard coded
TRGB: from table or from Mtrgb
Av, Dmod: From match
RGB region: Color limits by hand, TRGB / completeness or offset
Write optical LF
write ratio table
"""
from .pop_synth.stellar_pops import limiting_mag, rgb_agb_regions
import ResolvedStellarPops as rsp
import numpy as np
import sys
import os

def initialize_inputs():
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
            'matchphot': None,
            'Mtrgb': -4,
            'comp_frac': 0.9,
            'object_sfr_file': None,
            'object_mass': None,
            'outfile_loc': os.getcwd(),
            'file_imf': None,
            'object_cutoffmass': None}


def main(inps, angst=False):
    #import pdb; pdb.set_trace()
    rsp.fileio.ensure_dir(inps.outfile_loc)
    inps.target, inps.filter1, inps.filter2 = rsp.parse_pipeline(inps.matchphot)
    if inps.offset is None:
        print('finding completeness fraction')
        inps.comp_mag1, inps.comp_mag2 = limiting_mag(inps.fake_file,
                                                      inps.comp_frac)

    if angst:
        print('code a call to angst tables')
    else:
        assert inps.sfh_file is not None, 'need Av, angst target, or sfh_file'
        msfh = rsp.match.utils.MatchSFH(inps.sfh_file)
        inps.Av = msfh.Av
        # Trilegal needs dmod0 not dmodV!
        inps.dmod = msfh.dmod #- msfh.Av
        if inps.trgb is None:
            inps.trgb = rsp.astronomy_utils.Mag2mag(inps.Mtrgb, inps.filter2,
                                                    inps.photsys,
                                                    dmod=inps.dmod, Av=inps.Av)

    mag1, mag2 = np.loadtxt(inps.matchphot, unpack=True)

    if inps.mag_faint is None:
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
    inps.object_sfr_file = os.path.join(inps.outfile_loc,
                                        '%s.sfr' % inps.target.lower())
    rsp.match.utils.process_match_sfh(inps.sfh_file, outfile=inps.object_sfr_file)

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
    gal_inp.write_params(inps.galaxy_input, rsp.trilegal.utils.galaxy_input_fmt())

    return inps


if __name__ == "__main__":
    inp_obj = rsp.fileio.InputParameters(default_dict=initialize_inputs())
    inp_obj.add_params(rsp.fileio.load_input(sys.argv[1]))
    inp_obj = main(inp_obj, angst=False)
    inp_obj.write_params(sys.argv[1])
