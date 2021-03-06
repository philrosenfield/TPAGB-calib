
import galaxy_tests
import os
import ResolvedStellarPops as rsp
from TPAGBparams import snap_src
import numpy as np
import logging
logger = logging.getLogger()
#from sfh_tests_multi_proc import load_default_ancient_galaxies

def parse_sfh_data(filename, file_origin, frac=0.2):
    '''
    parse match sfh into a np.recarray
    ARGS:
    filename: a match sfh file to parse, needs to have at least:
              dtype = [('lagei', '<f8'),
                       ('lagef', '<f8'),
                       ('sfr', '<f8'),
                       ('sfr_errp', '<f8'),
                       ('sfr_errm', '<f8'),
                       ('mh', '<f8'),
                       ('mh_errp', '<f8'),
                       ('mh_errm', '<f8'),
                       ('mh_disp', '<f8')]
    file_origin: 'match' or 'match-grid'
        if 'match': reads the binned sfh as takes sfr_errs as is
        if 'match-grid': will overwrite sfr_errs.
            if sfr = 0: sfr_errp = min(sfr) * frac
            if sfr != 0: sfr_errp = sfr_errm = sfr * frac
    frac: multiplicitve factor to set sfr_err, only used if file_origin
        is set to match-grid.
    RETURNS:
    np.recarray of the sfh file

    use file_origin='match' only if Hybrid MonteCarlo has been run.
    '''
    if 'match-grid' == file_origin.lower() or 'match-hmc' == file_origin.lower():
        data = rsp.match_utils.read_binned_sfh(filename)
    elif 'match-old' == file_origin.lower():
        data = rsp.match_utils.read_match_old(filename)
    else:
        logger.error('please add a new data reader')

    if 'grid' in file_origin.lower():
        print 'CONSTANT %f SFR ERROR' % frac
        # at least have a uniform error that is 20% the smallest sfr.
        data.sfr_errp = np.min(data.sfr[data.sfr > 0]) * frac
        # now take 10% as nominal error in each sfr bin.
        data.sfr_errp[data.sfr > 0] = data.sfr[data.sfr > 0] * frac
        data.sfr_errm[data.sfr > 0] = data.sfr[data.sfr > 0] * frac
    return data


def completeness_table(targets, comp_val, ast_kw=None, outfile='default'):
    '''
    ex: completeness_table('ancients', 0.9)
    '''
    targets = galaxy_tests.load_targets(targets)
    if outfile == 'default':
        outfile = os.path.join(snap_src + '/tables/',
                               'completeness_%.2f.dat' % comp_val)

    fmt = '%(target)s %(opt_filter1)s %(opt_filter2)s %(ir_filter1)s %(ir_filter2)s \n'
    with open(outfile, 'w') as out:
        out.write('# completeness fraction: %.2f \n' % comp_val)
        out.write('# '\
                  + fmt.replace('%(', '').replace(')s', '').replace('\'', ''))
        for target in targets:
            print target
            ast_dict = find_completeness(target, comp_val)
            ast_dict['target'] = target
            out.write(fmt % ast_dict)
    logger.info('wrote %s' % outfile)


def find_completeness(target, comp_val, ast_kw=None):
    ast_kw = ast_kw or {}
    ast_kw = dict({'combined_filters': True,
                   'interpolate': True}.items() + ast_kw.items())
    ast_dict = {}

    fake_files = galaxy_tests.get_fake_files(target)
    for band, fake_file in zip(['opt', 'ir'], fake_files):
        asts = rsp.Galaxies.artificial_star_tests(fake_file)
        asts.completeness(**ast_kw)
        (ast_dict['%s_filter1' % band], ast_dict['%s_filter2' % band]) = \
            asts.get_completeness_fraction(comp_val)
    return ast_dict

def completeness_below_trgb():
    targets = galaxy_tests.ancients()
    ast_kw = {'combined_filters': True, 'interpolate': True}

    for target in targets:
        fake_files = galaxy_tests.get_fake_files(target)
        for band, fake_file in zip(['opt', 'ir'], fake_files):
            if band == 'opt':
                fits_src = snap_src + '/data/angst_no_trim'
            else:
                continue
                fits_src = 'default'
            gal = galaxy_tests.load_galaxy(target, band=band,
                                           fits_src=fits_src)
            gal.make_hess(binsize=0.1, hess_kw={'cbinsize': 0.05})
            
            asts = rsp.Galaxies.artificial_star_tests(fake_file)
            asts.completeness(**ast_kw)
            
            imag, icol = np.argwhere(gal.hess[2] == gal.hess[2].max())[0]
            rc_mag2 = gal.hess[1][imag]
            comp = asts.fcomp2(rc_mag2 - 1.)
            print '%s %s %.3f %.3f' % (target, band, rc_mag2, comp)

def completeness_corrections(dmag=0.1):
    '''
    get the completeness fraction for a given list of magnitudes.
    mag_bins can be either mag1 or mag2, will return both completeness
    fractions.
    '''
    ags = load_default_ancient_galaxies()
    ast_kw = {'combined_filters': True, 'interpolate': True}
    ast_dict = {}

    for i, target in enumerate(ags.data.target):
        fake_files = galaxy_tests.get_fake_files(target)
        ast_dict[target] = {}
        for band, fake_file in zip(['opt', 'ir'], fake_files):
            mag_bins = np.arange(ags.data['%s_min' % band][i],
                                 ags.data['%s_max' % band][i], dmag)
            asts = rsp.Galaxies.artificial_star_tests(fake_file)
            asts.completeness(**ast_kw)
            ast_c = asts.fcomp2(mag_bins)
            ast_dict[target]['%s_correction' % band] = ast_c
            ast_dict[target]['%s_bins' % band] = mag_bins

    return ast_dict

def write_completeness_corrections(outfile='default'):
    '''
    each line has target [band]_bins vals or target [band]_correction vals...
    '''
    if outfile == 'default':
        outfile = snap_src + '/tables/completeness_corrections.dat'
    ast_dict = completeness_corrections()
    with open(outfile, 'w') as out:
        for target, data_dict  in ast_dict.items():
            for key, vals in data_dict.items():
                out.write('%s %s ' % (target, key))
                out.write('%s \n' % ' '.join(['%g' % v for v in vals]))
    return

def read_completeness_corrections(filename='default'):
    '''
    format is
    target string floats
    like
    ddo71 opt_bin ....
    '''
    if filename == 'default':
        filename = snap_src + '/tables/completeness_corrections.dat'
    with open(filename, 'r') as f:
        lines = f.readlines()
    ast_dict = {}
    for line in lines:
        lsplit = line.strip().split()
        if not lsplit[0] in ast_dict.keys():
            ast_dict[lsplit[0]] = {}
        ast_dict[lsplit[0]][lsplit[1]] = np.array(lsplit[2:], dtype=float)
    return ast_dict

def read_completeness_table(table='default', absmag=False, uncertainties=False):
    if table == 'default':
        table = snap_src + '/tables/completeness_0.90.dat'

    if absmag is True:
        table = table.replace('.dat', '_absmag.dat')

    if uncertainties is True:
        table = table.replace('.dat', '_uncertainties.dat')
        dtype = [('target', '|S16'), ('opt_filter1', '<f8'),
                 ('opt_filter2', '<f8'), ('opt_color', '<f8'),
                 ('ir_filter1', '<f8'), ('ir_filter2', '<f8'),
                 ('ir_color', '<f8'), ]
    else:
        dtype = [('target', '|S16'), ('opt_filter1', '<f8'),
                 ('opt_filter2', '<f8'), ('ir_filter1', '<f8'),
                 ('ir_filter2', '<f8')]

    data = np.genfromtxt(table, dtype=dtype)
    return data
