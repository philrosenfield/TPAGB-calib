from optparse import OptionParser
import operator
import brewer2mpl
import ResolvedStellarPops as rsp
import fileIO
import LFUtils
import mk_sims
from multiprocessing import Pool
import itertools
import os
import sys
import numpy as np
#import pdb; pdb.set_trace()
import matplotlib.pyplot as plt
import GoogleSitesTable as gst
import matplotlib.nxutils as nxutils
import logging
logger = logging.getLogger()


# default locations.
tpcalib_dir = '/Users/phil/research/TP-AGBcalib/'
snap_src = os.path.join(tpcalib_dir, 'SNAP')

table_src = os.path.join(snap_src, 'tables')
plt_dir = os.path.join(snap_src, 'plots')
model_src = os.path.join(snap_src, 'models')
data_src = os.path.join(snap_src, 'data')
output_src = os.path.join(model_src, 'ast')
fits_src = os.path.join(data_src, 'galaxies')
fake_dir = os.path.join(data_src, 'fakes')
sfr_dir = os.path.join(data_src, 'sfh')
angst_data = rsp.angst_tables.AngstTables()


def edit_sfh(ID, new_sfr_dir, orig_sfr_dir=None, **changes):
    '''
    add or subtract a value from each bin in sfr file.

    possible changes:
    age
    sfr
    z
    feh (dex)
    zdisp - - not yet.

    example:
    IDs = all_IDs()
    new_sfr_dir = os.path.join(data_src, 'sfh - 0.3dex')
    new_sfrs = [edit_sfh(ID, new_sfr_dir, **{'FeH': - 0.3}) for ID in IDs]
    '''
    z_sun = 0.01524
    sfr_file = mk_sims.get_sfrFILE(ID, sfr_dir=orig_sfr_dir)
    orig_sfr = np.genfromtxt(sfr_file, names=['age', 'sfr', 'z'])
    new_sfr = orig_sfr.copy()
    for key in changes.keys():
        if key in new_sfr.dtype.names:
            new_sfr[key] += changes[key]
        elif key.lower() == 'feh':
            z = new_sfr['z']
            feh = np.log10(z / z_sun)
            new_feh = feh + changes[key]
            new_sfr['z'] = z_sun * 10 ** new_feh
        else:
            print '%s not found in %s' % (key, sfr_file)
            return -1

    rsp.fileIO.ensure_dir(new_sfr_dir)
    new_sfr_file = os.path.join(new_sfr_dir, os.path.split(sfr_file)[1])
    np.savetxt(new_sfr_file, new_sfr, fmt='%2.6f')
    return new_sfr_file


def tag_cmds(IDs):
    '''
    only need to do this once, hopefully. Just went though and added integer
    tags to each fits table so that it would be like trilegal output with the
    - l option.
    '''
    seqs = ['MS', 'RHeB', 'RGB']
    if type(IDs) == str:
        search_string = '*' + '*'.join((IDs, 'IR', '.fits'))
        fits = rsp.fileIO.get_files(fits_src, search_string)
    else:
        search_strings = ['*' + '*'.join((ID, 'IR', '.fits')) for ID in IDs]
        fits = np.squeeze([rsp.fileIO.get_files(fits_src, ss)
                           for ss in search_strings])
    rsp.annotate_cmd.define_color_mag_region(fits, seqs)
    return


def run_all(IDs, models):
    pool = Pool()

    res = []
    for ID, model in itertools.product(IDs, models):
        res.append(pool.apply_async(main, (ID, model),))

    for r in res:
        r.get()

    return


def get_trgb_fitsname(ID, band):
    if band == 'opt':
        fits = rsp.fileIO.get_files(fits_src,
                                    '*' + '*'.join((ID, 'trim', '.fits')))[0]
        (filter1, filter2) = os.path.split(fits)[1].split('.')[0].split('_')[-2:]
        trgb = angst_data.get_tab5_trgb_av_dmod(ID, filter1, filter2)[0]
    elif band == 'ir':
        # read data
        fits, = rsp.fileIO.get_files(fits_src,
                                     '*' + '*'.join((ID, 'IR', '.fits')))
        trgb = LFUtils.get_trgb_ir_nAGB(ID)[0]
    else:
        print 'photsys not found... check get_trgb_fitsname'
        sys.exit()
    return trgb, fits


def load_galaxy(ID, band):
    '''
    '''
    (trgb, fitsname) = get_trgb_fitsname(ID, band)

    kwargs = {'hla': False, 'trgb': trgb, 'photsys': 'wfc3snap', 'angst': True}
    kwargs['z'] = LFUtils.get_key_fromtable(ID, 'Z')

    try:
        taggedname = rsp.fileIO.replace_ext(fitsname, '.dat')
        gal = rsp.Galaxies.galaxy(taggedname, **kwargs)
    except:
        print 'no tagged data for %s' % fitsname
        gal = rsp.Galaxies.galaxy(fitsname, filetype='fitstable', **kwargs)

    return gal


def get_mix_modelname(model):
    mix = model.replace('cmd_input', '').split('.')[0].split('_')[0]
    model_name = '_'.join(model.split('.')[0].split('_')[1:])
    return mix, model_name


def get_fake_file(ID):
    return rsp.fileIO.get_files(fake_dir, '*%s*' %
                                ID.replace('C-0', 'C-').replace('C-', 'C'))[0]


def get_sfr_file(ID):
    return rsp.fileIO.get_files(sfr_dir, '*%s*' % ID)[0]


def match_metallicities(IDs):
    sfhs = [mk_sims.get_sfrFILE(ID) for ID in IDs]
    zs = {}
    for ID, sfh in zip(IDs, sfhs):
        age, sfr, z = np.loadtxt(sfh, unpack=True)
        zs[ID] = {'z': z, 'avez': np.mean(z)}
    return zs


def compare_metallicities():
    IDs = all_IDs()
    zs = match_metallicities(IDs)
    for ID in IDs:
        zs[ID]['zmeas'] = LFUtils.get_key_fromtable(ID, 'Z')

    for id, zdict in sorted(zs.items(), key=lambda(k, v): (v['avez'], k)):
        print '%s %.4f %.4f' % (id, zdict['avez'], zdict['zmeas'])


def get_stage_inds(fits, stage_name):
    inds, = np.nonzero(fits['stage'] == rsp.TrilegalUtils.get_stage_label(stage_name))
    return inds


def setup_trilegal(gal, model, object_mass=5e9):
    agb_model = model.replace('cmd_input_', '').replace('.dat', '').lower()

    table = '/Users/phil/research/TP-AGBcalib/SNAP/tables/table.dat'
    with open(table, 'r') as t:
        lines = t.readlines()
    for line in lines:
        line = line.strip().split()
        if line[0] == gal.target:
            mag_lim_filt = line[4]
            mag_lim = float(line[5])

    object_dist = 10 ** ((5 + gal.dmod) / 5)
    sfr_file = get_sfr_file(gal.target)
    file_mag = 'tab_mag_odfnew/tab_mag_%s.dat' % gal.photsys
    galinp_kw = {'object_sfr_file': get_sfr_file(gal.target),
                 'file_mag': file_mag,
                 'mag_limit_val': mag_lim,
                 'mag_num': rsp.TrilegalUtils.find_mag_num(file_mag, mag_lim_filt),
                 'object_sfr_file': sfr_file,
                 'object_av': gal.Av,
                 'object_dist': object_dist,
                 'object_mass': object_mass,
                 'object_sfr_mult_factorA': 1e9}

    galaxy_input = os.path.join(snap_src, 'input', 'input_%s.dat' % gal.target)
    rsp.TrilegalUtils.change_galaxy_input(galaxy_input, **galinp_kw)
    trilegal_output = os.path.join(snap_src, 'output', 'output_%s_%s.dat' %
                                   (gal.target, agb_model))

    return galaxy_input, trilegal_output, galinp_kw


def setup_data_normalization(gal, maglims):
    # find normalization region
    # (when hand picking cmd_regions, RGB includes AGB!)
    gal.irgb = gal.stage_inds('RGB')
    # the stars in the data marked as rgb between mag lims.
    rgb_norm = rsp.math_utils.between(gal.mag2, maglims[0], maglims[1],
                                      inds=gal.irgb)
    points = np.column_stack((gal.color[rgb_norm], gal.mag2[rgb_norm]))
    # the 1 std around color mean
    cmean = np.mean(gal.color[rgb_norm])
    cstd = np.std(gal.color[rgb_norm])
    verts = np.array([[cmean - cstd, maglims[0]],
                      [cmean - cstd, maglims[1]],
                      [cmean + cstd, maglims[1]],
                      [cmean + cstd, maglims[0]],
                      [cmean - cstd, maglims[0]]])

    ndata_stars = np.nonzero(nxutils.points_inside_poly(points, verts))[0].size
    return verts, ndata_stars


def make_normalized_simulation(ID, model, photsys='wfc3snap', over_write=False,
                               object_mass=5e6, maglims='trgb', band='ir',
                               offsets=(1.5, 0.5), mk_sims_args={},
                               count_offset=0.):
    '''
    Wrapper for running trilegal and LFUtils.calc_LF
    Will continue to run trilegal until galaxy is high enough mass to have
    proper normalization.

    many attributes of galaxy and simgalaxy instance are set in
    LFUtils.calc_LF.

    inputs:
    ID of target: string
    cmd_input model to run trilegal: string

    optional params:
    xxxover_write: bool [False]
    xxx    force running of trilegal (will run if no output file exists)
    maglims: list, tuple or string. ['trgb']
        either the dim, bright maglims or 'trgb'
    offsets: list, tuple: (dim, bright) default: [(1.5, 0.5)]
        must exist if maglims='trgb'

        maglims will be trgb + dim, trgb + bright
    object_mass: float [5e6]
        initial mass to run trilegal.

    returns:
        gal instance
        sgal instance
        ks test p_value
        maglims
    '''
    gal = load_galaxy(ID, band)
    galaxy_input, trilegal_output, galinp_kw = setup_trilegal(gal, model,
                                                              object_mass=object_mass)
    cmd_input = os.path.join('/Users/phil/research/padova_apps/cmd_inputfiles/', model)
    fake_file = get_fake_file(ID)

    # set maglim if not already set.
    if maglims == 'trgb':
        maglims = (gal.trgb + offsets[0], gal.trgb + offsets[1])

    verts, ndata_stars = setup_data_normalization(gal, maglims)

    #gal.iagb = rsp.math_utils.brighter(gal.mag2, gal.trgb, inds=gal.irgb)

    # initializations
    go = 0
    normalization = 1e9
    # arbitrary normalization factor set to 0.75... 1.0 is upper limit,
    # because that would mean there are the same number of stars in the
    # normalization region for both the sim and the data. Any number lower than
    # one is just a bit of overkill so you can randomly draw and get a good
    # stat. sample.
    while normalization > .75:
        # object_mass increase factor is hard coded, could be kwarg.
        if go > 0:
            object_mass *= 5.
            galinp_kw['object_mass'] = object_mass
            rsp.TrilegalUtils.change_galaxy_input(galaxy_input, **galinp_kw)
        go += 1

        # run trilegal
        print 'Trying %s %s, Mass=%g, Attempt %i' % (ID, model, object_mass, go)

        rsp.TrilegalUtils.run_trilegal(cmd_input, galaxy_input, trilegal_output)

        # load sim galaxy
        sgal = rsp.Galaxies.simgalaxy(trilegal_output, gal.filter1, gal.filter2,
                                      count_offset=count_offset)
        # "correct" for asts
        rsp.Galaxies.ast_correct_trilegal_sim(sgal, fake_file=fake_file)
        sgal.load_ast_corrections()

        sgal.normalize('rgb', useasts=True, ndata_stars=ndata_stars, by_stage=False,
                       verts=verts)
        # calcuate the LFs
        #p_value = LFUtils.calc_LF(gal, sgal, maglims)
        #print p_value

        # update normalization
        normalization = sgal.rgb_norm
        print 'normalization', normalization

    sgal.ID = ID
    sgal.model = model
    sgal.mix_modelname(model)

    if over_write:
        sgal.object_mass = object_mass
        print '\n %s necessary object_mass = %g\n' % (gal.name, object_mass)

    return


def all_IDs():
    IDs = ["SCL-DE1",
           "NGC2403-HALO-6",
           "NGC7793-HALO-6",
           "UGC-04459",
           "UGC-4305-1",
           "UGC-4305-2",
           "UGC-5139",
           "DDO78",
           "DDO82",
           "KDG73",
           "KKH37",
           "NGC0300-WIDE1",
           "NGC404",
           "NGC2403-DEEP",
           "NGC2976-DEEP",
           "NGC3741",
           "NGC4163",
           "UGC8508",
           "UGCA292",
           "NGC3077-PHOENIX",
           "IC2574-SGS",
           "HS117",
           "DDO71"]
    return IDs


def main(IDs, models, band='ir', mk_sims_args={}, make_plots=False, publish_plots=False,
         count_offset=0., outfile='result_tab.dat'):
    '''
    big wrapper for galaxy tests

    kwargs:
    make_plots: LFUtils.plotLFIR and .diagnostic_cmd
    publish_plots: call GoogleSitesTable
    '''

    if os.path.isfile(outfile):
        print outfile, 'exists. appending.'
        out = open(outfile, 'a')
    else:
        out = open(outfile, 'w')
        out.write('# ID model p_value NRGB_data NAGB_data NRGB_model NAGB_model\n')

    for ID, model in itertools.product(IDs, models):
        sim_kwargs = {'mk_sims_args': mk_sims_args, 'band': band}
        gal, sgal, p_value, maglims = make_normalized_simulation(ID, model,
                                                                 count_offset=count_offset,
                                                                 **sim_kwargs)
        #write_spread_catalog(sgal, **mk_sims_args)

        if make_plots:
            fig_loc = os.path.join(plt_dir, sgal.mix, sgal.model_name)
            diag_loc = os.path.join(fig_loc, 'diag')
            rsp.fileIO.ensure_dir(diag_loc + '/')

            figname = os.path.join(diag_loc,
                                   '%s_%s_diag.png' % (ID, sgal.model_name))
            LFUtils.diagnostic_cmd(sgal, gal.trgb, figname=figname,
                                   inds=sgal.rel_ind)
            # should take fig_loc...
            LFUtils.plot_LFIR(gal, sgal, p_value, maglims)

        out.write('%s %s %.3f %i %i %i %i\n' % (sgal.ID,
                                                sgal.model,
                                                p_value,
                                                gal.irgb.size,
                                                gal.iagb.size,
                                                sgal.rel_rgb.size,
                                                sgal.rel_agb.size))

    out.close()

    for model in models:
        #compare_sims(IDs, model)
        if publish_plots:
            html_file = os.path.join(diag_loc, sgal.model_name + '_diag.html')
            gst.side_by_side(diag_loc, sgal.model_name, html_file)

            html_file = os.path.join(fig_loc, sgal.model_name + '.html')
            gst.one_col(fig_loc, sgal.model_name, html_file)

    return


def load_galaxies(IDs, model, **kwargs):
    sgals = []
    gals = []
    #IDs = all_IDs()
    out_dir = kwargs.get('out_dir', model_src)
    #for ID, model in itertools.product(IDs, models):
    for ID in IDs:
        spread_file = os.path.join(out_dir, 'spread',
                                   'spread_output_%s_model_%s' % (ID, model))
        gal = rsp.Galaxy.galaxy(ID, 'ir')
        sgal = rsp.Galaxy.simgalaxy(spread_file, gal.filter1, gal.filter2)
        sgal.ID = ID
        sgal.model = model
        sgal.mix_modelname(model)
        sgal.get_fits()
        #maglims = (np.nan, np.nan)
        #p_value = LFUtils.calc_LF(gal, sgal, maglims, normalize=False)
        sgal.file_ast = mk_sims.get_fakFILE(ID)
        sgals.append(sgal)
        gals.append(gal)

    # galaxies - ify
    SGals = rsp.Galaxy.galaxies(sgals)
    Gals = rsp.Galaxy.galaxies(gals)
    return Gals, SGals


def load_galaxies_by_z(IDs, model):
    Gals, SGals = load_galaxies(IDs, model)
    # select on the ones with measured z:
    galsz = Gals.finite_key('z')
    # sort by z.
    galsz_sort = sorted(galsz, key=lambda galaxy: galaxy.z)

    # use that to grab the sim gals
    sgalsz = []
    for g in galsz_sort:
        sgalsz += [sg for sg in SGals.galaxies if sg.target == g.target]

    zs = [g.z for g in galsz_sort]
    return zs, galsz_sort, sgalsz


def compare_sims(IDs, model, fig_name=None):
    zs, galsz_sort, sgalsz = load_galaxies_by_z(IDs, model)
    ax = plt.axes()
    for z, g, s in zip(zs, galsz_sort, sgalsz):
        ax.plot(z, float(s.rel_agb.size) / float(g.iagb.size),
                'o', ms=5, color='black')
        bright_rgb = LFUtils.brighter(s.ast_mag2, g.trgb, inds=s.irgb)
        ax.plot(z,
                float((s.rel_agb.size - bright_rgb.size)) / float(g.iagb.size),
                'o', ms=5, color='red')
        ax.annotate(g.target,
                    xy=(z, float((s.rel_agb.size - bright_rgb.size)) /
                        float(g.iagb.size) + 0.01), ha='center')
        ax.annotate(g.target,
                    xy=(z, float(s.rel_agb.size) / float(g.iagb.size) + 0.01),
                    ha='center')
    ax.set_ylabel(r'$N_{AGB, model} / N_{AGB, data}$', fontsize=20)
    ax.set_xlabel(r'$Z$', fontsize=20)
    ax.set_title(r'$\rm{%s}$' % s.model_name.replace('_', '\ '))
    rsp.fileIO.ensure_dir(os.path.join(plt_dir, s.mix))
    if not fig_name:
        fig_name = os.path.join(plt_dir, s.mix, s.model_name + '.png')
    plt.savefig(fig_name)
    plt.close()
    print 'wrote %s' % fig_name
    return


def chi2_plot(IDs, models):
    ax = plt.axes()
    cols = brewer2mpl.get_map('Dark2', 'qualitative', len(models)).mpl_colors
    l = 0.
    for i, model in enumerate(models):
        l += .004
        zs, galsz_sort, sgalsz = load_galaxies_by_z(IDs, model)
        for z, g, s in zip(zs, galsz_sort, sgalsz):
            ax.plot(z, s.chi2, 'o', ms=5, color=cols[i])
            ax.annotate(g.target, xy=(z + 0.0001, s.chi2 + 0.01), ha='left',
                        va='center', fontsize='8')
        ax.annotate('$%s$' % s.model_name.replace('_', '\ '), xy=(l, 3),
                    color=cols[i], fontsize=20)
    ax.set_ylabel(r'$\chi^2$', fontsize=20)
    ax.set_xlabel(r'$Z$', fontsize=20)
    fig_name = os.path.join(plt_dir, s.mix, 'chi2.png')
    ax.set_ylim(2, 18)
    plt.savefig(fig_name)
    plt.close()
    return


def write_spread_catalog(sgal, outfile=None, **kwargs):
    '''
    writes a slice of the trilegal output catalog.
    slice is on ast recovered stars that are randomly selected by
    normalization.
    i.e sgal.rec and sgal.rel_ind.
    returns string outfile name.
    '''
    out_dir = kwargs.get('outdir')
    if out_dir is None:
        out_dir = model_src
    if outfile is None:
        outfile = os.path.join(out_dir, 'spread', sgal.name.replace('ast',
                                                                    'spread'))
        rsp.fileIO.ensure_dir(outfile)
    if not os.path.isfile(outfile):
        out = open(outfile, 'w')
        # get the header
        sort_keys = sorted(sgal.data.key_dict.iteritems(),
                           key=operator.itemgetter(1))
        col_keys = [k[0] for k in sort_keys]

        # slice data on recovered inds and rel_inds.
        dslice = np.squeeze(sgal.data.data_array[[sgal.rec], :])[sgal.rel_ind]

        out.write('# {}\n'.format(' '.join(col_keys)))
        np.savetxt(out, dslice, fmt='%-7.6f')
        print 'wrote %s' % outfile
    else:
        print '%s exists. Write spread catalog will not overwrite.' % outfile

    return outfile


def setup_match(name_fmt, bg_name=None, phot_name=None, **kwargs):
    '''
    sets up the match directory structure
    '''
    if bg_name is not None and phot_name is not None:
        mdirs = {'msgs': '.msg',
                 'output': '.fit',
                 'plots': '',
                 'pars': '.par',
                 'bg': bg_name,
                 'phot': phot_name}
        match_dir = os.path.join(model_src, 'match')
    else:
        match_dir = kwargs.get('match_dir')
        kwargs.pop('match_dir')
        mdirs = kwargs

    mdirdict = {}
    for d, ext in mdirs.items():
        new_dir = os.path.join(match_dir, d + '/')
        rsp.fileIO.ensure_dir(new_dir)
        # if no extension, return the plot directory
        if len(ext) == 0:
            mdirdict[d] = new_dir
            continue
        # if extension less than 5 chars, use name_fmt.ext
        if len(ext) <= 5:
            mdirdict[d] = os.path.join(new_dir, '%s%s' % (name_fmt, ext))
        # if more than 5 characters, assume it's a file name.
        else:
            mdirdict[d] = os.path.join(new_dir, ext)

    return mdirdict


def match_tests(IDs, model):
    Gals, SGals = load_galaxies(IDs, model)
    for g, s in zip(Gals.galaxies, SGals.galaxies):
        name_fmt = '%s_%s_%s' % (s.ID, s.mix, s.model_name)
        bg_name = s.name.replace('spread', 'bg')
        phot_name = rsp.fileIO.replace_ext(os.path.split(g.name)[1], '.match')
        match_dict = setup_match(name_fmt, bg_name, phot_name)
        match_dict['plots'] = os.path.join(match_dict['plots'], s.mix,
                                           s.model_name)
        rsp.fileIO.ensure_dir(match_dict['plots'] + '/')

        match_bg = rsp.match_utils.make_match_bg_cmd(s.ast_mag1, s.ast_mag2,
                                                     outfile=match_dict['bg'])
        print match_bg
        #phot = rsp.match_utils.make_match_bg_cmd(g.mag1, g.mag2,
        #                                    outfile=match_dict['phot'])

        match_kwargs = {'dmod': s.data.get_col('m - M0')[0],
                        'Av': s.data.get_col('Av')[0],
                        'filter1': s.filter1.replace('F', 'IR'),
                        'filter2': s.filter2.replace('F', 'IR'),
                        'color': s.ast_color,
                        'mag': s.ast_mag2,
                        'pmfile': match_dict['pars']}

        pm_file = rsp.match_utils.make_calcsfh_param_file(match_bg, **match_kwargs)
        # run match

        match_out = rsp.match_utils.call_match(pm_file, match_dict['phot'],
                                               s.file_ast, match_dict['output'],
                                               match_dict['msgs'])
        chi, fit = rsp.match_utils.get_fit(match_out)
        # make plot
        alabel = '$\mathrm{%s}' % s.model_name.replace('_', '\ ')
        metallicity = g.z
        if not np.isfinite(g.z):
            metallicity = '...'
            alabel += '\ Z = %s$' % metallicity
        else:
            alabel += '\ Z = %.4f$' % metallicity
        cmdgrid = match_out + '.cmd'
        labs = ['$%s$' % s.ID, alabel, '$\mathrm{Diff}$',
                '$\chi^2=%g\ \mathrm{Sig}$' % chi]
        grid = rsp.match_utils.pgcmd(cmdgrid, labels=labs)
        grid[2].set_ylabel(r'$%s$' % g.filter2, fontsize=20)
        grid[2].set_xlabel(r'$%s - %s$' % (g.filter1, g.filter2), fontsize=20)
        figname = os.path.join(match_dict['plots'], os.path.split(cmdgrid)[1])
        print figname
        plt.savefig('%s.png' % figname)
        plt.close()
        s.chi = chi
        s.fit = fit
    return Gals, SGals


def short_list():
    IDs = ['DDO82', 'NGC2403-HALO-6', 'NGC2976-DEEP', 'NGC4163',
           'NGC7793-HALO-6', 'UGC8508', 'UGCA292']
    return IDs


if __name__ == "__main__":

    ch = logging.StreamHandler()
    ch.setLevel(logging.DEBUG)
    logger.addHandler(ch)
    logger.info('start of run')

    parser = OptionParser()

    usage = "%prog model [options]"

    parser = OptionParser(usage=usage)

    parser.add_option("-s", action="store_true", default=False,
                      help="run on short set of galaxies")

    parser.add_option("-m", action="store_true", default=False,
                      help="make plots")

    parser.add_option("-p", action="store_true", default=False,
                      help="publish plots")

    parser.add_option("-o", action="store_true", default=False,
                      help="out directory set to SNAP/models")

    parser.add_option('-i', action='store_true', default=False,
                      help='providing AGBTracksUtils inputfile')

    (options, args) = parser.parse_args()

    if options.s:
        IDs = short_list()[0]
        models = [args[0]]

    kwargs = {}
    if not options.i:
        if options.m:
            kwargs['make_plots'] = True
        if options.p:
            kwargs['publish_plots'] = True
        if options.o:
            kwargs['outdir'] = '/Users/phil/research/TP-AGBcalib/SNAP/'
            kwargs['count_offset'] = 1.
    else:
        print 'reading directions from input file'
        input_file = args[0]
        infile = fileIO.input_file(input_file)
        IDs = infile.IDs
        if IDs is None:
            'No galaxies listed, doing short list.'
            IDs = short_list()
        agb_mix = infile.agb_mix
        set_name = infile.set_name
        track_set = '_'.join((agb_mix, set_name))
        models = ['%s.dat' % track_set]
        kwargs['outdir'] = infile.galaxy_outdir
        if infile.google_table:
            kwargs['make_plots'] = True
            kwargs['publish_plots'] = True

    main(IDs=IDs, models=models, **kwargs)

    '''
    #IDs = [IDs[0]]
    #model = "cmd_input_CAF09_S_SCS.dat"

    models = ["cmd_input_CAF09_s_aug12.dat",
              "cmd_input_CAF09_S_SCS.dat",
              "cmd_input_CAF09_S_SCSFG.dat",
              "cmd_input_CAF09_S_SCSFG_ETA2.dat"]
    #for model in models:
    #    match_tests(IDs, model)
    #    compare_sims(IDs, model)
    chi2_plot(IDs, models)
    sfr_dir = os.path.join(data_src, 'sfh - 0.3dex')
    out_dir = os.path.join(model_src, '0.3dex')
    #plt_dir ...
    kwargs = {'make_plots': True,
              'publish_plots': True,
              'sfr_dir': sfr_dir,
              'out_dir': out_dir}
    main(IDs=IDs, **kwargs)
    '''

    '''
    for model in models:
        sgal, gal = compare_sims(IDs, model)
        fig_loc = os.path.join(plt_dir, sgal.mix, sgal.model_name)
        diag_loc = os.path.join(fig_loc, 'diag')
        html_file = os.path.join(diag_loc, sgal.model_name + '_diag.html')
        gst.side_by_side(diag_loc, sgal.model_name, html_file)
        html_file = os.path.join(fig_loc, sgal.model_name + '.html')
        gst.one_col(fig_loc, sgal.model_name, html_file)
    '''
