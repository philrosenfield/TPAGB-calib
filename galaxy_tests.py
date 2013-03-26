from optparse import OptionParser
import operator
import brewer2mpl
import ResolvedStellarPops as rsp
import fileIO
import LFUtils
import mk_sims
import os
import numpy as np
#import pdb; pdb.set_trace()
import matplotlib.pyplot as plt
import matplotlib.nxutils as nxutils
import logging
logger = logging.getLogger()


# default locations.
snap_src = '/Users/phil/research/TP-AGBcalib/SNAP'

plt_dir = os.path.join(snap_src, 'plots')
model_src = os.path.join(snap_src, 'models')

angst_data = rsp.angst_tables.AngstTables()


def tag_cmds(IDs):
    '''
    only need to do this once, hopefully. Just went though and added integer
    tags to each fits table so that it would be like trilegal output with the
    - l option.
    '''
    seqs = ['MS', 'RHeB', 'RGB']
    fits_src = os.path.join(snap_src, 'data', 'galaxies')
    if type(IDs) == str:
        search_string = '*' + '*'.join((IDs, 'IR', '.fits'))
        fits = rsp.fileIO.get_files(fits_src, search_string)
    else:
        search_strings = ['*' + '*'.join((ID, 'IR', '.fits')) for ID in IDs]
        fits = np.squeeze([rsp.fileIO.get_files(fits_src, ss)
                           for ss in search_strings])
    rsp.annotate_cmd.define_color_mag_region(fits, seqs)
    return


def load_galaxy(ID, band='ir'):
    '''
    '''
    fits_src = os.path.join(snap_src, 'data', 'galaxies')
    if band is None:
        fits_src = os.path.join(snap_src, 'data', 'opt_ir_matched_v2')
    fitsname = rsp.fileIO.get_files(fits_src, '*%s*fits' % ID)
    if len(fitsname) > 1:
        filetype = 'fitstable'
    else:
        filetype = 'agbsnap'
    fitsname = ir_or_opt_file(fitsname, band=band)
    gal_kw = {'hla': False, 'photsys': 'wfc3snap', 'angst': True,
              'filetype': filetype, 'band': band}

    gal_kw['z'] = LFUtils.get_key_fromtable(ID, 'Z')

    gal = rsp.Galaxies.galaxy(fitsname, **gal_kw)

    return gal


def ir_or_opt_file(filenames, band='ir'):
    '''
    if filenames is a list len 1: does nothing.
    if band is ir returns the file with IR in the filename,
    else returns the other file, or the first file without IR in it.
    '''
    if len(filenames) <= 1:
        filename, = filenames
    else:
        file_ir = [a for a in filenames if 'IR' in a][0]
        file_opt = [a for a in filenames if not 'IR' in a][0]
        if band == 'ir':
            filename = file_ir
        else:
            filename = file_opt
    return filename


def get_mix_modelname(model):
    mix = model.replace('cmd_input', '').split('.')[0].split('_')[0]
    model_name = '_'.join(model.split('.')[0].split('_')[1:])
    return mix, model_name


def get_fake_files(ID, band=None):
    fake_dir = os.path.join(snap_src, 'data', 'fakes')
    fake_files = rsp.fileIO.get_files(fake_dir, '*%s*' %
                                      ID.replace('C-0', 'C-').replace('C-', 'C'))
    file_ir = [a for a in fake_files if 'IR' in a][0]
    file_opt = [a for a in fake_files if not 'IR' in a][0]
    if band is None:
        fake_file = [file_opt, file_ir]
    elif band == 'opt':
        fake_file = file_opt
    else:
        fake_file = file_ir
    return fake_file


def get_sfr_file(ID):
    sfr_dir = os.path.join(snap_src, 'data', 'sfh')
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


def setup_trilegal(gal, model, object_mass=5e9):
    '''
    Sets up files for trilegal simulations (same files that will be overwritten
        in a loop).

    1) Writes an initial galaxy_input file ()
    2) Creates the kwargs to pass to
    rsp.TrilegalUtils.change_galaxy_input to change the galaxy_input file
    in the loop.
    3) Creates the string name of the simulation output file (trilegal_output)
    Sets the 50% completeness mag as the mag limit for the trilegal simulation

    Possible source of future error:
    The sfr_file that are used as TPAGB input have time bins in Gyr.
    Change object_sfr_mult_factorA if that statement is no longer true.

    globals
    snap_src file locations

    input
    gal, galaxy object with attributes:
    gal.photsys
    gal.dmod
    gal.target
    gal.Av
    model, agb model. string, e.g: 'cmd_input_CAF09_S_SCS.dat'

    '''
    agb_model = model.replace('cmd_input_', '').replace('.dat', '').lower()

    object_dist = 10 ** ((5 + gal.dmod) / 5)
    sfr_file = get_sfr_file(gal.target)
    file_mag = 'tab_mag_odfnew/tab_mag_%s.dat' % gal.photsys

    galinp_kw = {'object_sfr_file': get_sfr_file(gal.target),
                 'file_mag': file_mag,
                 'mag_limit_val': gal.comp50mag2 + 3.,
                 'mag_num': rsp.TrilegalUtils.find_mag_num(file_mag,
                                                           gal.filter2),
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


def setup_data_normalization(gal, filt1, filt2, band='ir', leo_method=False):
    '''
    returns vertices and number of data stars used in normalization.

    prepares for rsp.Galaxies.simgalaxy.normalize by creating vertices of a
    cmd-polygon to use for normalization and then finding the number of data
    stars in that polygon.

    polygon is a box centered at the mean color of rgb stars from the data.
    The rgb region has been hand picked for the galaxy object, it must have
    been a "taggedfits" file.
    the color edges of the polygon are +- 1/2 std from mean color.
    the mag edges of the polygon are the maglims, which are typically set to
    some mag below the trgb.
    '''
    assert hasattr(gal, 'maglims'), \
        'maglims must be an attribute of galaxy object.'

    f1 = gal.filters.index(filt1) + 1
    f2 = gal.filters.index(filt2) + 1
    mag1 = gal.__getattribute__('mag%i' % f1)
    mag2 = gal.__getattribute__('mag%i' % f2)
    color = mag1 - mag2

    if leo_method is False:
        # find normalization region
        # (when hand picking cmd_regions, RGB includes AGB!)
        gal.irgb = gal.stage_inds('RGB')
        # the stars in the data marked as rgb between mag lims.
        rgb_norm = rsp.math_utils.between(mag2, gal.maglims[0],
                                          gal.maglims[1],
                                          inds=gal.irgb)
        # the 1/2 std around color mean
        cmean = np.mean(color[rgb_norm])
        cstd = np.std(color[rgb_norm])/2.
        col_min = cmean - cstd
        col_max = cmean + cstd
    else:
        # Leo's method is to use all the stars between the mag lims.

        rgb_norm = rsp.math_utils.between(mag2, gal.maglims[0], gal.maglims[1])
        col_min = np.min(color)
        col_max = np.max(color)

    verts = np.array([[col_min, gal.maglims[0]],
                      [col_min, gal.maglims[1]],
                      [col_max, gal.maglims[1]],
                      [col_max, gal.maglims[0]],
                      [col_min, gal.maglims[0]]])

    points = np.column_stack((color, mag2))
    inds, = np.nonzero(nxutils.points_inside_poly(points, verts))
    gal.rgb_norm_inds = list(set(rgb_norm) & set(inds))
    ndata_stars = len(gal.rgb_norm_inds)

    return verts, ndata_stars


def make_normalized_simulation(ID, model, filt1, filt2, photsys='wfc3snap',
                               over_write=False, object_mass=5e6,
                               maglims='trgb', band='ir', offsets=(1.5, 0.5),
                               count_offset=0., mass_inc_fact=5,
                               run_trilegal=True, norm_threshold=0.75,
                               leo_method=False, spread_outfile=None):
    '''
    Will continue to run trilegal until galaxy is high enough mass to have
    proper normalization.

    Arbitrary normalization factor set to 0.75... 1.0 is upper limit,
    because that would mean there are the same number of stars in the
    normalization region for both the sim and the data. Any number lower than
    one is just a bit of overkill so you can randomly draw and get a good
    stat. sample.

    inputs:
    ID of target: string
    cmd_input model to run trilegal: string

    optional params:

    band [ir] ir or opt switch, which data file to be read in.

    run_trilegal [True]
        do the simulation (else just load the existing one)

    maglims: float list or tuple else a string. ['trgb']
        either the dim, bright maglims or 'trgb'
        if trgb, uses offsets.

    offsets: list, tuple: (dim, bright) default: [(1.5, 0.5)]
        must exist if maglims='trgb'

        maglims will be trgb + dim, trgb + bright

    object_mass: float [5e6]
        initial mass to run trilegal.

    returns:
        gal instance
            with more attributes!
            gal.irgb
            gal.rgb_norm_inds
        sgal instance with ast corrections
            with more attributes!
            sgal.rgb_norm
            sgal.rgb_norm_inds
            sgal.object_mass
            sgal.norm_verts
            sgal.target
            sgal.model
    '''
    gal = load_galaxy(ID, band=band)
    trgb = gal.trgb
    if filt1 == 'default':
        filt1 = gal.filter1
        filt2 = gal.filter2
    if band == 'ir' and hasattr(gal, 'mag3'):
        trgb = gal.ir_trgb

    galaxy_input, trilegal_output, galinp_kw = setup_trilegal(gal, model,
                                                              object_mass=object_mass)

    cmd_input = os.path.join('/Users/phil/research/padova_apps/cmd_inputfiles/', model)
    fake_files = get_fake_files(ID, band=band)

    # set maglim if not already set.
    if maglims == 'trgb':
        gal.maglims = (trgb + offsets[0], trgb + offsets[1])

    verts, ndata_stars = setup_data_normalization(gal, filt1, filt2,
                                                  leo_method=leo_method)

    # initializations
    go = 0
    normalization = 1e9

    while normalization > norm_threshold:
        if go > 0:
            # if we've been through this already, increase mass and try again.
            object_mass *= mass_inc_fact
            galinp_kw['object_mass'] = object_mass
            rsp.TrilegalUtils.change_galaxy_input(galaxy_input, **galinp_kw)
        go += 1

        # run trilegal
        logger.debug('Trying %s %s, Mass=%g, Attempt %i' %
                     (ID, model, object_mass, go))

        if run_trilegal is True:
            rsp.TrilegalUtils.run_trilegal(cmd_input, galaxy_input,
                                           trilegal_output)

        # load sim galaxy
        sgal = rsp.Galaxies.simgalaxy(trilegal_output, gal.filter1, gal.filter2,
                                      count_offset=count_offset)
        # "correct" for asts
        if leo_method is True:
            # Leo's method will write to a file.
            # There must be filters in the file name for spread_angst_ir to
            # work properly.
            if spread_outfile is None:
                spread_outfile = os.path.join(sgal.base, 'ast_%s_%s_%s' %
                                              (gal.filter1, gal.filter2, sgal.name))
            spread_outfile2 = None
            if hasattr('gal', 'filter3'):
                spread_outfile2 = spread_outfile.replace(sgal.name, '%s_%s_%s' %
                                                        (gal.filter3, gal.filter4,
                                                         sgal.name))

        rsp.Galaxies.ast_correct_trilegal_sim(sgal, fake_file=fake_files,
                                              leo_method=leo_method,
                                              spread_outfile=spread_outfile,
                                              spread_outfile2=spread_outfile2)
        if spread_outfile2 is not None:
            spread_outfile = spread_outfile2

        if leo_method is True:
            # now load the ast file we just wrote
            # TODO: open up spread_angst and send the array back!
            sgal = rsp.Galaxies.simgalaxy(spread_outfile, gal.filter1, gal.filter2,
                                          count_offset=count_offset)

        ast_check = sgal.load_ast_corrections()
        assert ast_check == 1, 'problem with ast corrections'

        sgal.normalize('rgb', filt1, filt2, useasts=True, by_stage=False,
                       ndata_stars=ndata_stars, verts=verts)

        # update normalization
        normalization = sgal.rgb_norm
        print 'normalization', normalization

    sgal.norm_inds = sgal.rgb_norm_inds
    sgal.target = ID
    sgal.model = model
    sgal.mix_modelname(model)
    sgal.object_mass = object_mass
    sgal.norm_verts = verts
    return sgal, gal


def gi10_overlap():
    return ["DDO78", "DDO71", "SCL-DE1"]


def all_IDs():
    return ["SCL-DE1",
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


def repeat_girardi10(maglims='trgb'):
    filt1 = 'default'
    filt2 = 'default'
    band = 'opt'
    targets = ['DDO78', 'DDO71', 'SCL-DE1']

    models = ['cmd_input_CAF09_S_JAN13.dat'
              'cmd_input_gi10_rev.dat',
              'cmd_input_gi10_rev_old_tracks.dat']

    offsets = (2., 0.)
    norm_sim_kw = {'offsets': offsets, 'band': band, 'leo_method': True,
                   'object_mass': 1e8}

    plot_LF_kw = {'ylim': (28, 20), 'xlim2': (0.8, 4e4)}

    for model in models:
        if model == 'cmd_input_gi10_rev.dat':
            model_title = 'Gi10\ S12ND'
        elif model == 'cmd_input_gi10_rev_old_tracks.dat':
            model_title = 'Gi10'
        elif model == 'cmd_input_CAF09_S_JAN13.dat':
            model_title = 'JAN13'

        plot_LF_kw['model_title'] = model_title

        for target in targets:
            if target == 'DDO78':
                norm_sim_kw['object_mass'] = 5e8

            sgal, gal = make_normalized_simulation(target, model, filt1, filt2,
                                                   **norm_sim_kw)
            if gal.filter1 == 'F475W':
                plot_LF_kw['xlim'] = (-0.75, 4)
            else:
                plot_LF_kw['xlim'] = (-0.75, 3)

            smg = rsp.Galaxies.sim_and_gal(gal, sgal)
            smg.nrgb_nagb(band=band)
            smg.make_LF(filt1, filt2, plot_LF_kw=plot_LF_kw, comp50=True)


def ir_rgb_agb_ratio():
    filt1 = 'default'
    filt2 = 'default'
    band = 'ir'
    targets = ['DDO78', 'DDO71', 'SCL-DE1']
    models = ['cmd_input_CAF09_S_JAN13.dat',
              'cmd_input_gi10_rev.dat',
              'cmd_input_gi10_rev_old_tracks.dat']

    offsets = (1.5, 0.)
    norm_sim_kw = {'offsets': offsets, 'band': band, 'leo_method': True,
                   'object_mass': 1e8, 'run_trilegal': False}

    plot_LF_kw = {'ylim': (25.5, 20), 'xlim': (1, 2.5), 'xlim2': (0.8, 4e4)}

    for model in models:
        if model == 'cmd_input_gi10_rev.dat':
            model_title = 'Gi10\ S12ND'
        elif model == 'cmd_input_gi10_rev_old_tracks.dat':
            model_title = 'Gi10'
        elif model == 'cmd_input_CAF09_S_JAN13.dat':
            model_title = 'JAN13'

        plot_LF_kw['model_title'] = model_title

        for target in targets:
            if target == 'DDO78':
                norm_sim_kw['object_mass'] = 5e8

            sgal, gal = make_normalized_simulation(target, model, filt1, filt2,
                                                   **norm_sim_kw)

            smg = rsp.Galaxies.sim_and_gal(gal, sgal)
            smg.nrgb_nagb(band=band)
            smg.make_LF(gal.filter1, gal.filter2, plot_LF_kw=plot_LF_kw,
                        comp50=True)


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

    #main(IDs=IDs, models=models, **kwargs)
