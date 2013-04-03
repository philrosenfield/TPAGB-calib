import logging
logging.basicConfig(filename='galaxy_tests.log',level=logging.DEBUG)
logger = logging.getLogger()
logger.info('start of run')
import ResolvedStellarPops as rsp
import LFUtils
import os
import sys
import difflib
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.nxutils as nxutils

angst_data = rsp.angst_tables.AngstTables()


def sort_angst_data_table(table, targets):
    '''
    order table in the same index order as target list.
    table  (np.array) must have key 'target' which matches values of targets
    '''        
    # targets as known by snap table
    snap_targets = [difflib.get_close_matches(t, table['target'])[0]
                    for t in targets]
    
    snap_tab_sort = [list(table['target']).index(st) for st in snap_targets]
    return table[snap_tab_sort]

    
def gaussian_mixture_modeling():
    '''
    This is silly.
    '''
    from scipy.stats import norm
    from sklearn.mixture import GMM
    from astroML.datasets import fetch_sdss_sspp
    from astroML.decorators import pickle_results
    from astroML.plotting.tools import draw_ellipse

    targets = all_targets()
    # load galaxies class
    gals = rsp.Galaxies.galaxies([load_galaxy(t, band='ir') for t in targets])
    # combine all galaxies
    gals.squish('Color', 'Mag2', 'Trgb')
    mean_trgb = np.mean(gals.Trgbs)
    offset = gals.Trgbs - np.mean(gals.Trgbs)
    gals.Mag2o = np.concatenate([g.Mag2 - offset[i] for i, g in enumerate(gals.galaxies)])
    X = np.vstack([gals.Colors, gals.Mag2o]).T

    N = np.arange(1, 14)

    @pickle_results("GMM_CMD.pkl")
    def compute_GMM(N, covariance_type='full', n_iter=1000):
        models = [None for n in N]
        for i in range(len(N)):
            print N[i]
            models[i] = GMM(n_components=N[i], n_iter=n_iter,
                            covariance_type=covariance_type)
            models[i].fit(X)
        return models

    models = compute_GMM(N)
    AIC = [m.aic(X) for m in models]
    BIC = [m.bic(X) for m in models]
    i_best = np.argmin(BIC)
    gmm_best = models[i_best]
    print "best fit converged:", gmm_best.converged_
    print "BIC: n_components =  %i" % N[i_best]
    fig, ax = plt.subplots()
    colorbins = 501
    magbins = 501
    H, colorbins, magbins = np.histogram2d(gals.Colors, gals.Mag2o, (colorbins, magbins))
    Xgrid = np.array(map(np.ravel, np.meshgrid(0.5*(colorbins[:-1]+colorbins[1:]), 0.5*(magbins[:-1]+magbins[1:])))).T
    log_dens = gmm_best.score(Xgrid).reshape((501, 501))
    fig, ax = plt.subplots()
    ax.imshow(np.exp(log_dens),
              origin='lower', interpolation='nearest', aspect='auto',
              extent=[colorbins[0], colorbins[-1],
                      magbins[-1], magbins[0]], cmap=plt.cm.binary)

    ax.scatter(gmm_best.means_[:, 0], gmm_best.means_[:, 1], c='r')
    for mu, C, w in zip(gmm_best.means_, gmm_best.covars_, gmm_best.weights_):
        draw_ellipse(mu, C, scales=[1.5], ax=ax, fc='none', ec='k')
    

def multi_galaxy_hess():
    targets = all_targets()
    # load galaxies class
    gals = rsp.Galaxies.galaxies([load_galaxy(t, band='ir') for t in targets])
    # combine all galaxies
    gals.squish('Color', 'Mag2', 'Trgb')

    mean_trgb = np.mean(gals.Trgbs)
    offset = gals.Trgbs - np.mean(gals.Trgbs)
    gals.Mag2o = np.concatenate([g.Mag2 - offset[i] for i, g in enumerate(gals.galaxies)])
    hess_kw = {'binsize': 0.1, 'cbinsize': 0.05}
    hess = rsp.astronomy_utils.hess(gals.Colors, gals.Mag2o, **hess_kw)
    # now split by mean color.
    # snap table 3 in the same order of the targets list
    snap_tab3 = sort_angst_data_table(angst_data.snap_tab3, targets)

    mean_color = snap_tab3['mean_color']
    bluer, = np.nonzero(mean_color < 0.861)
    redder, = np.nonzero(mean_color >= 0.861)

    gals.squish('Color', **{'inds': bluer, 'new_attrs': ['bColor']})

    gals.squish('Color', **{'inds': redder, 'new_attrs': ['rColor']})

    gals.bMag2o = np.concatenate([g.Mag2 - offset[i]
                                  for i, g in enumerate(gals.galaxies[bluer])])

    gals.rMag2o = np.concatenate([g.Mag2 - offset[i]
                                  for i, g in enumerate(gals.galaxies[redder])])

    
    blue_hess = rsp.astronomy_utils.hess(gals.bColor, gals.bMag2o, **hess_kw)
    red_hess = rsp.astronomy_utils.hess(gals.rColor, gals.rMag2o, **hess_kw)

    # offset = Trgbs-np.mean(Trgbs)
    # Mag2o = np.concatenate([g.Mag2 - offset[i] for i, g in enumerate(gals)])
    # fig, ax = gals[0].plot_cmd(Color, Mag2o, threshold=100, levels=5)
    # fig.colorbar(gals[0].cs)


def agb_rheb_separation():
    '''
    right now this is done by hand. :)
    '''
    targets = all_targets()
    #targets.pop(targets.index('NGC404'))
    gals = [load_galaxy(t, band='ir') for t in targets]
    # Color = np.concatenate([g.Color for g in gals])
    # Mag1 = np.concatenate([g.Mag1 for g in gals])
    # Mag2 = np.concatenate([g.Mag2 for g in gals])
    Trgbs = np.array([g.Trgb for g in gals])
    mean_trgb = np.mean(Trgbs)
    # offset = Trgbs-np.mean(Trgbs)
    # Mag2o = np.concatenate([g.Mag2 - offset[i] for i, g in enumerate(gals)])
    # fig, ax = gals[0].plot_cmd(Color, Mag2o, threshold=100, levels=5)
    # fig.colorbar(gals[0].cs)
    # here's the place where there will be a double gaussian fitting.
    verts = np.array([[0.7, mean_trgb],
                      [0.77, -6.6],
                      [0.83, -7.4],
                      [0.83, -10],
                      [5, -10],
                      [5, mean_trgb],
                      [0.7, mean_trgb]])
    # so it's going to be this gal[i] agb's are within: verts[:,1] + offset[i]
    # [g.name, offset[i] for i, g in enumerate(gals)]
    return verts


def load_agb_verts(gal, leo_method=False):
    '''
    there is a file that contains the contents of this dictionary, it's vert_file
    I just didn't feel like writing a reader for it so I pasted it. What?
    don't look at me like that I have a lot to do.
    vert_file = '/Users/phil/research/TP-AGBcalib/code/TPAGB-calib/agb_rheb_sep.dat'
    '''
    if leo_method is False:
        offset_dict = {'SCL-DE1': 0.186,
                       'NGC2403-HALO-6': -0.155,
                       'NGC7793-HALO-6': -0.307,
                       'UGC-04459': 0.063,
                       'UGC-4305-1': 0.002,
                       'UGC-4305-2': 0.012,
                       'UGC-5139': 0.105,
                       'DDO78': -0.054,
                       'DDO82': -0.103,
                       'KDG73': 0.157,
                       'KKH37': 0.085,
                       'NGC0300-WIDE1': -0.212,
                       'NGC2403-DEEP': -0.089,
                       'NGC2976-DEEP': -0.220,
                       'NGC3741': 0.143,
                       'NGC4163': 0.003,
                       'UGC8508': 0.077,
                       'UGCA292': 0.394,
                       'NGC3077-PHOENIX': -0.297,
                       'IC2574-SGS': 0.044,
                       'HS117': 0.092,
                       'DDO71': 0.075}
        Verts = agb_rheb_separation()
        verts = Verts.copy()
        verts[:, 1] += offset_dict[gal.target]

        vColor = verts[:, 0]
        vMag2 = verts[:, 1]
        vMag1 = vColor + vMag2
        Mag2mag_kw = {'Av': gal.Av, 'dmod': gal.dmod}
        vmag1 = rsp.astronomy_utils.Mag2mag(vMag1, gal.filter1, gal.photsys,
                                            **Mag2mag_kw)

        vmag2 = rsp.astronomy_utils.Mag2mag(vMag2, gal.filter2, gal.photsys,
                                            **Mag2mag_kw)

        verts = np.column_stack((vmag1 - vmag2, vmag2))
    else:
        magdim = gal.trgb
        magbright = 0.
        colmin = np.min(gal.color)
        colmax = np.max(gal.color)
        verts = np.array([[colmin, magdim],
                          [colmin, magbright],
                          [colmax, magbright],
                          [colmax, magdim],
                          [colmin, magdim]])
    return verts


def load_sim_masses(ID):
    '''
    I let the code loop overnight and took these masses that were necessary for
    jan13 tracks to run. Some are a bit higher than necessary, but by no more
    than a factor of 2.5
    '''
    mass_dict = {'SCL-DE1': 1e+08,
                 'NGC2403-HALO-6': 1e+08,
                 'NGC7793-HALO-6': 1e+08,
                 'UGC-04459': 1e+08,
                 'UGC-4305-1': 5e+08,
                 'UGC-4305-2': 5e+08,
                 'UGC-5139': 1e+08,
                 'DDO78': 5e+08,
                 'DDO82': 2.5e+09,
                 'KDG73': 2e+07,
                 'KKH37': 1e+08,
                 'NGC0300-WIDE1': 1e+08,
                 'NGC404': 5e+08,
                 'NGC2403-DEEP': 1e+08,
                 'NGC2976-DEEP': 5e+08,
                 'NGC3741': 1e+08,
                 'NGC4163': 5e+08,
                 'UGC8508': 1e+08,
                 'UGCA292': 2e+07,
                 'NGC3077-PHOENIX': 1e+08,
                 'IC2574-SGS': 2.5e+09,
                 'HS117': 2e+07,
                 'DDO71': 1e+08}
    return mass_dict[ID]


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
        try:
            file_opt = [a for a in filenames if not 'IR' in a][0]
        except IndexError:
            print 'looks like no optical fits table'
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
    try:
        file_opt = [a for a in fake_files if not 'IR' in a][0]
    except IndexError:
        print 'looks like no optical fake file for %s.' % ID
        file_opt = None
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


def compare_metallicities():
    IDs = all_targets()
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


def setup_data_normalization(gal, filt1, filt2, band=None, leo_method=False):
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
        if gal.irgb is None:
            gal.irgb, = np.nonzero(gal.color > 0.4)
            logger.info('normalization color: using the stdev around the mean of stars in mag range using color cut of 0.4')
        else:
            logger.info('normalization color: using the stdev around the mean of stars of tagged rgb stars in mag range')

        # the stars in the data marked as rgb between mag lims.
        rgb_norm = rsp.math_utils.between(mag2, gal.maglims[0], gal.maglims[1],
                                          inds=gal.irgb)
        # the 1/2 std around color mean
        cmean = np.mean(color[rgb_norm])
        cstd = np.std(color[rgb_norm])/2.
        col_min = cmean - cstd
        col_max = cmean + cstd
    else:
        # Leo's method is to use all the stars between the mag lims.
        logger.info('normalization color: using all stars in mag range')
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
    gal.norm_verts = verts
    logger.debug('norm verts:')
    logger.debug(verts)
    logger.info('data stars in normalization range: %i' % ndata_stars)
    return verts, ndata_stars


def make_normalized_simulation(gal, model, filt1, filt2, photsys='wfc3snap',
                               over_write=False, object_mass=5e6,
                               maglims='trgb', band='ir', offsets=(1.5, 0.5),
                               count_offset=0., mass_inc_fact=5,
                               run_trilegal=True, norm_threshold=0.75,
                               leo_method=False, spread_outfile=None,
                               leo_norm=False, leo_ast=False,
                               spread_outfile2=None, trilegal_output=None,
                               norm_fname=None):
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

    if leo_method is True:
        leo_norm = True
        leo_ast = True

    if filt1 == 'gal':
        filt1 = gal.filter1
        filt2 = gal.filter2
        logger.info('filter1 set to %s' % filt1)
        logger.info('filter2 set to %s' % filt2)

    # only in the case of 4-filter will there be an attr ir_trgb.
    # otherwise it's gal.trgb which is always of gal.filter2.
    if band == 'ir' and hasattr(gal, 'mag3'):
        trgb = gal.ir_trgb
        logger.info('using IR trgb')
    else:
        trgb = gal.trgb
        logger.info('using %s trgb' % gal.filter2)

    if run_trilegal is True or trilegal_output is None:
        galaxy_input, trilegal_output, galinp_kw = setup_trilegal(gal, model,
                                                                  object_mass=object_mass)

        cmd_input = os.path.join('/Users/phil/research/padova_apps/cmd_inputfiles/', model)

    # might not be needed, but better to load them outside the while loop.
    fake_files = get_fake_files(gal.target, band=band)

    # set maglim if not already set.
    if maglims == 'trgb':
        gal.maglims = (trgb + offsets[0], trgb + offsets[1])
        logger.info('normalization mag: %f below the trgb to %f above' % (offsets[0], offsets[1]))
    else:
        logger.info('normalization mag: %f to %f' % (offsets[0], offsets[1]))

    verts, ndata_stars = setup_data_normalization(gal, filt1, filt2,
                                                  leo_method=leo_norm)

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

        if run_trilegal is True:
            # run trilegal
            logger.debug('Trying %s %s, Mass=%g, Attempt %i' %
                         (gal.target, model, object_mass, go))

            rsp.TrilegalUtils.run_trilegal(cmd_input, galaxy_input,
                                           trilegal_output)

        # load sim galaxy
        sgal = rsp.Galaxies.simgalaxy(trilegal_output, gal.filter1, gal.filter2,
                                      count_offset=count_offset)
        ast_check = sgal.load_ast_corrections()
        if not hasattr(sgal, 'ast_mag2'):
            # "correct" for asts
            if leo_ast is True:
                # Leo's method will write to a file.
                if spread_outfile is None:
                    spread_outfile = os.path.join(sgal.base, 'ast_%s_%s_%s' %
                                                  (gal.filter1, gal.filter2, 
                                                   sgal.name))
                else:
                    assert gal.filter1 in spread_outfile, \
                    'There must be filters in the file name for spread_angst_ir to work properly.'

                if hasattr('gal', 'filter3'):
                    spread_outfile2 = spread_outfile.replace(sgal.name,
                                                             '%s_%s_%s' %
                                                            (gal.filter3,
                                                             gal.filter4,
                                                             sgal.name))

            rsp.Galaxies.ast_correct_trilegal_sim(sgal, fake_file=fake_files,
                                                  leo_method=leo_ast,
                                                  spread_outfile=spread_outfile,
                                                  spread_outfile2=spread_outfile2)
            if spread_outfile2 is not None:
                spread_outfile = spread_outfile2

            if leo_ast is True:
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
        logger.debug('normalization %.4f' % normalization)

    if norm_fname is None:
        norm_fname = os.path.join(sgal.base, sgal.name.replace('.dat','_inds.dat'))
    write_norm_inds(norm_fname, sgal)
    sgal.norm_inds = sgal.rgb_norm_inds
    sgal.target = gal.target
    sgal.model = model
    sgal.mix_modelname(model)
    sgal.object_mass = object_mass
    sgal.norm_verts = verts
    return sgal


def write_norm_inds(filename, sgal):
    with open(filename, 'a') as f:
        f.write('# %s:\n' % os.path.join(sgal.base, sgal.name))
        f.write('%s \n' % ', '.join(map(str, sgal.rgb_norm_inds)))
    logger.info('normalization indices written to %s' % filename)


def load_normalized_simulation(target, model, band=None, input_file=None,
                               maglims=None, offsets=None, leo_norm=False):

    gal = load_galaxy(target, band=band)
    sim_file = load_ast_file(gal, model)
    if input_file is not None:
        inputs = rsp.fileIO.load_inputs(input_file)
        offsets = inputs['offsets']
        maglims = inputs['maglims']
        leo_norm = inputs['leo_norm']

    gal.maglims = maglim_offsets(offsets, maglims, gal.trgb)

    verts, ndata_stars = setup_data_normalization(gal, gal.filter1, gal.filter2,
                                                  leo_method=leo_norm)
    
    # just takes the first... it's a quick load.
    sgal = rsp.Galaxies.simgalaxy(sim_file, filter1=gal.filter1, filter2=gal.filter2)
    sgal.load_ast_corrections()

    sgal.target = gal.target
    sgal.model = model
    sgal.mix_modelname(model)
    sgal.norm_verts = verts
    return gal, sgal

def call_mc_norm(do_all=True, gi10=False):
    if do_all is True:
        targets = all_targets()
        models = ['cmd_input_CAF09_S_JAN13.dat',
                  'cmd_input_gi10_rev_old_tracks.dat',
                  'cmd_input_gi10_rev.dat']
        mc_norm_kw = {'band': 'ir', 'offsets': (1.5, 0.), 'leo_norm': False,
                      'maglims': 'trgb'}
    
    if gi10 is True:
        targets = gi10_overlap()
        models = ['cmd_input_gi10_rev_old_tracks.dat']
        mc_norm_kw = {'band': 'opt', 'offsets': (2., 0.), 'leo_norm': True,
                      'maglims': 'trgb', 'leo_method': True}
    for target in targets:
        for model in models:
            mc_norm(target, model, **mc_norm_kw)
            #sgal_rgb_agb(target, model, **mc_norm_kw)

def mc_norm(target, model, band=None, input_file=None, maglims=None,
            offsets=None, leo_norm=False, leo_method=False):
    
    gal, sgal = load_normalized_simulation(target, model, band=band, 
                                           input_file=input_file,
                                           maglims=maglims, offsets=offsets, 
                                           leo_norm=leo_norm)

    rgb_verts, ndata_rgb = setup_data_normalization(gal, gal.filter1, gal.filter2,
                                                    leo_method=leo_norm)

    points = np.column_stack((gal.color, gal.mag2))
    spoints = np.column_stack((sgal.ast_color[sgal.rec], sgal.ast_mag2[sgal.rec]))

    reg_inds, = np.nonzero(nxutils.points_inside_poly(spoints, rgb_verts))
    nsim_stars = len(reg_inds)

    agb_verts = load_agb_verts(gal, leo_method=leo_method)
    ndata_agb = len(np.nonzero(nxutils.points_inside_poly(points, agb_verts))[0])

    nrgb_nagb_data = float(ndata_agb) / float(ndata_rgb)
    nmodel = np.array([], dtype=float)

    normalization = float(ndata_rgb) / float(nsim_stars)
    for i in range(10001):
        # random sample the data distribution
        rands = np.random.random(len(sgal.rec))
        ind, = np.nonzero(rands < normalization)
        
        # the number of sim stars in the rgb box set by data verts
        srgb_norm, = np.nonzero(nxutils.points_inside_poly(spoints[ind],
                                                           rgb_verts))
        nsim_rgb = float(len(srgb_norm))

        sagb_norm, =  np.nonzero(nxutils.points_inside_poly(spoints[ind],
                                                            agb_verts))
        nsim_agb = float(len(sagb_norm))

        nrgb_nagb_sim = nsim_agb / nsim_rgb
        nmodel = np.append(nmodel, nrgb_nagb_sim)
    print '%s %.3f %.3f %.4f %s' % (target, nrgb_nagb_data, np.mean(nmodel),
                                    np.std(nmodel), model)

def sgal_rgb_agb(target, model, band=None, input_file=None, maglims=None,
                 offsets=None, leo_norm=False):

    gal, sgal = load_normalized_simulation(target, model, band=band, 
                                           input_file=input_file,
                                           maglims=maglims, offsets=offsets, 
                                           leo_norm=leo_norm)
    sgal.all_stages('RGB', 'TPAGB')
    rgb_inds = list(set(sgal.irgb) & set(sgal.rec))
    agb_inds = list(set(sgal.itpagb) & set(sgal.rec))

    raw_ratio = float(sgal.itpagb.size)/float(sgal.irgb.size)
    rgb_verts = sgal.norm_verts
    agb_verts = load_agb_verts(gal)

    spoints = np.column_stack((sgal.ast_color, sgal.ast_mag2))

    srgb_norm, = np.nonzero(nxutils.points_inside_poly(spoints[rgb_inds],
                                                       rgb_verts))
    nsim_rgb = float(len(srgb_norm))

    sagb_norm, =  np.nonzero(nxutils.points_inside_poly(spoints[agb_inds],
                                                        agb_verts))
    nsim_agb = float(len(sagb_norm))    
    nrgb_nagb_sim = nsim_agb / nsim_rgb
    
    # again with no asts.
    sapoints = np.column_stack((sgal.color, sgal.mag2))

    sargb_norm, = np.nonzero(nxutils.points_inside_poly(sapoints[sgal.irgb],
                                                        rgb_verts))
    nasim_rgb = float(len(sargb_norm))
    saagb_norm, =  np.nonzero(nxutils.points_inside_poly(sapoints[sgal.itpagb],
                                                         agb_verts))
    nasim_agb = float(len(saagb_norm))    
    ratio_no_ast = nasim_agb / nasim_rgb
                                                       
    print '%s %.3f %.3f %.3f %s' % (target, nrgb_nagb_sim, raw_ratio, ratio_no_ast, model)

def compare_models():
    dtype = [('galaxy', 'S16'), ('data_ratio', '<f8'), ('JAN13', '<f8'), 
             ('JAN13_stdev', '<f8'),  ('Gi10', '<f8'), ('Gi10_stdev', '<f8'),
             ('Gi10_S12ND', '<f8'), ('Gi10_S12ND_stdev', '<f8'), ('Z', '<f8')]

    fyeah = np.genfromtxt('ir_rgb_agb_result_n10000.dat', dtype=dtype)

    models = ['JAN13', 'Gi10', 'Gi10_S12ND']
    resids = np.array([fyeah[m] - fyeah['data_ratio'] for m in models])
    from matplotlib.ticker import NullFormatter

    fig, ax = plt.subplots()
    [ax.plot([i,i,i], resids.T[i], color='black') for i in range(len(resids.T))]
    [ax.plot(resids[i], 'o', label='$%s$' % models[i].replace('_','\ ')) for i in range(len(models))]
    [ax.text(i, max(resids.T[i]+0.01), '$%s$' % fyeah['galaxy'][i], fontsize=10) for i in range(len(resids.T))]
    [ax.text(i, min(resids.T[i]-0.01), '$%.4f$' % fyeah['Z'][i], fontsize=10) for i in range(len(resids.T))]
    ax.xaxis.set_major_formatter(NullFormatter()) 

def gi10_overlap():
    return ['DDO78', 'DDO71', 'SCL-DE1']


def all_targets():
    return ['SCL-DE1',
            'DDO78',
            'DDO71',
            'NGC7793-HALO-6',
            'NGC0300-WIDE1',
            'NGC3077-PHOENIX',
            'NGC2403-DEEP',
            'NGC2403-HALO-6',
            'UGC-5139',
            'DDO82',
            'IC2574-SGS',
            'UGC-4305-1',
            'UGC-4305-2',
            'NGC4163',
            'UGC8508',
            'UGC-04459',
            'NGC3741',
            'UGCA292',
            'HS117',
            'KDG73',
            'KKH37',
            'NGC2976-DEEP']


def short_list():
    return['DDO82', 'NGC2403-HALO-6', 'NGC2976-DEEP', 'NGC4163',
           'NGC7793-HALO-6', 'UGC8508', 'UGCA292']


def load_ast_file(gal, model):
    model_short = model.replace('cmd_input_','').lower()
    search_str = '*'.join(('ast_%s' % gal.filter1, gal.filter2, gal.target,
                           model_short))

    sim_file, = rsp.fileIO.get_files(os.path.join(snap_src, 'output'),
                                     search_str)
    return sim_file

def maglim_offsets(offsets, maglims, trgb=None):
    if maglims == 'trgb':
        maglims = (trgb + offsets[0], trgb + offsets[1])
    return maglims


def read_norm_inds_file(sgal, filename=None):
    if filename is None:
        filename = os.path.join(sgal.base, sgal.name.replace('.dat', '_inds.dat'))

    with open(filename, 'r') as f:
        lines = [l for l in f.readlines() if not l.startswith('#')]

    norm_inds = []
    for line in lines:
        inds = np.array(map(int, line.strip().split(',')))
        norm_inds.append(inds)
    return norm_inds

def ir_rgb_agb_ratio(renormalize=False, filt1=None, filt2=None, band=None,
                     targets=None, run_trilegal=False, leo_ast=True, 
                     offsets=(1.5, 0.), models=None, maglims=None,
                     leo_method=False, leo_norm=False, make_plot=True,
                     xlim=None, ylim=None, xlim2=None, add_boxes=True):
    
    if targets == 'all':
        targets = all_targets()
    elif targets == 'gi10':
        targets = gi10_overlap()

    if type(targets) == str:
        targets = [targets]

    if type(models) == str:
        models = [models]

    norm_sim_kw = {'offsets': offsets, 'band': band, 'leo_ast': leo_ast,
                   'run_trilegal': run_trilegal, 'leo_method': leo_method,
                   'leo_norm': leo_norm}

    if make_plot is True:
        # these are chosen for IR plots.
        if xlim is None:
            xlim = (-0.5, 1.5)
        if ylim is None:
            ylim = (25.3, 18)
        if xlim2 is None:
            xlim2 = (0.8, 4e4)
        plot_LF_kw = {'ylim': ylim, 'xlim': xlim,
                      'xlim2': xlim2}

    result_dict = {}
    for target in targets:
        norm_sim_kw['object_mass'] = load_sim_masses(target)
        gal = load_galaxy(target, band=band)

        for model in models:
            logger.info('%s: %s' % (target, model))
            model_title = model.replace('.dat', '').split('_')[-1]
            if model == 'cmd_input_gi10_rev_old_tracks.dat':
                model_title = 'Gi10'
            elif model == 'cmd_input_gi10_rev.dat':
                model_title = 'Gi10\ S12ND'

            if renormalize is True:
                norm_sim_kw['trilegal_output'] = load_ast_file(gal, model)

            sgal = make_normalized_simulation(gal, model, filt1, filt2,
                                              **norm_sim_kw)
            agb_verts = load_agb_verts(gal, leo_method=leo_method)
            smg = rsp.Galaxies.sim_and_gal(gal, sgal)
            nrgb_nagb_data, nrgb_nagb_sim = smg.nrgb_nagb(band=band,
                                                          agb_verts=agb_verts)
            if make_plot is True:
                plot_LF_kw['model_title'] = model_title
                smg.make_LF(gal.filter1, gal.filter2, plot_LF_kw=plot_LF_kw,
                            comp50=True, add_boxes=add_boxes)
            
            sub_dict = {'data': nrgb_nagb_data, model_title: nrgb_nagb_sim}
            if result_dict.has_key(target):
                result_dict[target].update(sub_dict)
            else:
                result_dict[target] = sub_dict

    print '   ', ' '.join(result_dict[result_dict.keys()[0]].keys())
    for key, rd in result_dict.items():
        print key, ' '.join(['%.3f' % i for i in  rd.values()])

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


if __name__ == "__main__":
    #ch = logging.StreamHandler()
    #ch.setLevel(logging.DEBUG)
    #logger.addHandler(ch)

    #import cProfile
    #command = 'ir_rgb_agb_ratio()'
    #cProfile.runctx(command, globals(), locals())

    inputs = rsp.fileIO.load_input(sys.argv[1])

    if inputs['debug'] is True:
        import pdb
        pdb.set_trace()
        del inputs['debug']

    global snap_src
    snap_src = inputs['snap_src']

    del inputs['snap_src']

    ir_rgb_agb_ratio(**inputs)
else:
    snap_src = '/Users/phil/research/TP-AGBcalib/SNAP'
    
