import galaxy_tests
import os
import ResolvedStellarPops as rsp
from TPAGBparams import snap_src
import numpy as np
import logging
logger = logging.getLogger()


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
