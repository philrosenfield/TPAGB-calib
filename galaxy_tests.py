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
from matplotlib.ticker import MaxNLocator
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
from pprint import pprint
import brewer2mpl

angst_data = rsp.angst_tables.AngstTables()


def cmd_contamination(gal, dcolor=0.1, dmag=0.3, thresh=5):
    '''
    , verts_file_fmt=None,
    diag_plot=True, isoch_base=None,
    isoc_search_fmt=None, out_file_fmt=None,
    out_fig_fmt=None, M11=False, galaxy_table=None,
    rdcolor=None, rdmag=None
    '''
    # rough outline of the RGB/AGB area. Exclude MS (~0.0)
    color_limits = [0.4, 1.5]
    bright_mag_limit = np.min(gal.mag2)

    # minimum size of dmag bin spacing
    dmag = 0.5

    #inds, = np.nonzero((gal.mag2 < gal.trgb) &
    #                   (gal.color > color_limits[0]) &
    #                   (gal.color < color_limits[1]))
    #mag2 = gal.mag2[inds]
    #color = gal.color[inds]
    points = np.column_stack((gal.color, gal.mag2))
    mag_bins = np.arange(bright_mag_limit, gal.trgb + dmag,  dmag)

    for u, l in zip(mag_bins, np.roll(mag_bins,-1))[:-1]:
        verts = np.array([[color_limits[0], l],
                          [color_limits[0], u],
                          [color_limits[1], u],
                          [color_limits[1], l],
                          [color_limits[0], l]])
        inside, = np.nonzero(nxutils.points_inside_poly(points, verts))
        
        if len(inside) < thresh:
            print len(inside)
            continue

        rstars_in_agb, agb_in_rstars, poisson_noise, nstars, color_sep = \
            gal.double_gaussian_contamination(verts, diag_plot=True, dcol=0.1)
        lost_agb = agb_in_rstars
        contaminated_agb = rstars_in_agb
        print u, l, color_sep
        # still doesn't work... not picking up double gaussians....
def get_imf(target):
    filename = '/Users/phil/research/TP-AGBcalib/code/TPAGB-calib/best_fits_from_match_runs.dat'
    dtype = [('ID', '|S5'), ('Galaxy', '|S14'), ('filter1', '|S5'), ('filter2', '|S5'), ('Av', '<f8'), ('IMF', '<f8'), ('dmod', '<f8'), ('dlogZ', '<f8')]
    data = np.genfromtxt(filename, dtype=dtype)
    imf, = data['IMF'][np.nonzero(data['Galaxy']==target)]
    if imf == 1.30:
        tab_imf = 'tab_imf/imf_salpeter_match.dat'        
    elif imf == 1.35:
        tab_imf = 'tab_imf/imf_salpeter.dat'
    else:
        print 'imf problem.'
    return tab_imf


def make_sfh_plots():
    zctmps = ['/Users/phil/research/Italy/WFC3SNAP/PHIL/SFRfiles/noAGB/fullreszctmps/10915_DDO82_F606W_F814W.gst.sfh.zctmp',
              '/Users/phil/research/Italy/WFC3SNAP/PHIL/SFRfiles/noAGB/fullreszctmps/9755_IC2574-SGS_F555W_F814W.gst.match.sfh.zctmp',
              '/Users/phil/research/Italy/WFC3SNAP/PHIL/SFRfiles/noAGB/fullreszctmps/10605_UGC-4305-1_F555W_F814W.gst.sfh.zctmp',
              '/Users/phil/research/Italy/WFC3SNAP/PHIL/SFRfiles/noAGB/fullreszctmps/10605_UGC-4305-2_F555W_F814W.gst.sfh.zctmp',
              '/Users/phil/research/Italy/WFC3SNAP/PHIL/SFRfiles/noAGB/fullreszctmps/10915_NGC4163_F606W_F814W.gst.sfh.zctmp',
              '/Users/phil/research/Italy/WFC3SNAP/PHIL/SFRfiles/noAGB/fullreszctmps/10915_UGC8508_F475W_F814W.gst.match.sfh.zctmp'] 
    [rsp.match_utils.plot_zctmp(z) for z in zctmps]


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
    

def multi_galaxy_hess(targets=None, split_by_color=False, ax=None, imshow_kw={},
                      make_hess=False, band='ir'):
    if ax is None:
        fig, ax = plt.subplots()

    if targets is None:
        targets = all_targets()
    # load galaxies class
    gals = rsp.Galaxies.galaxies([load_galaxy(t, band=band) for t in targets])
    # combine all galaxies
    gals.squish('Color', 'Mag2', 'Trgb')

    if split_by_color is True:
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
        # not finished...
    else:
        mean_trgb = np.mean(gals.Trgbs)
        offset = gals.Trgbs - np.mean(gals.Trgbs)
        gals.Mag2o = np.concatenate([g.Mag2 - offset[i] for i, g in enumerate(gals.galaxies)])
        if make_hess is True:
            cmin, cmax, cbinsize = -0.25, 1.25, 0.01
            mmin, mmax, binsize = -10.5, -2, 0.05
            cbins = np.arange(cmin, cmax, cbinsize)
            mbins = np.arange(mmin, mmax, binsize)
            hess_kw = {'binsize': 0.1, 'cbinsize': 0.01, 'cbin': cbins, 'mbin': mbins}
            hess = rsp.astronomy_utils.hess(gals.Colors, gals.Mag2o, **hess_kw)
            imshow_kw = dict({'norm': LogNorm(vmin=None, vmax=hess[2].max()),
                              'cmap': plt.cm.gray_r, 'interpolation': 'nearest'}.items() + 
                              imshow_kw.items())
            ax = rsp.astronomy_utils.hess_plot(hess, imshow_kw=imshow_kw)
        else:
            fig, ax = gals.galaxies[0].plot_cmd(gals.Colors, gals.Mag2o,
                                                threshold=50, levels=5,
                                                hist_bin_res=.1)
            #fig.colorbar(gals.galaxies[0].cs)
        verts = agb_rheb_separation()
        ax.plot(verts[:, 0], verts[:, 1], lw=2, color='navy')
        ax.hlines(mean_trgb, *ax.get_xlim(), lw=2, color='red', zorder=100)
        ax.set_xlim(cmin, cmax)

    ax.set_xlabel('$%s-%s$' % (gals.filter1, gals.filter2), fontsize=20)
    ax.set_ylabel('$%s$' % gals.filter2, fontsize=20)
    ax.tick_params(labelsize=16)
    return fig, ax


def agb_rheb_separation(targets=None):
    '''
    right now this is done by hand. :)
    '''
    if targets is None:
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

    # I was a little off...
    verts[:, 0] += 0.04
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

    gal_dict_inp = {'photsys': gal.photsys,
                    'filter1': gal.filter2,
                    'object_mass': object_mass,
                    'object_sfr_file': get_sfr_file(gal.target)}

    gal_dict = rsp.TrilegalUtils.galaxy_input_dict(**gal_dict_inp)
    gal_inp = rsp.fileIO.input_parameters(default_dict=gal_dict)

    galinp_kw = {'mag_limit_val': gal.comp50mag2 + 3.,
                 'object_sfr_file': sfr_file,
                 'object_av': gal.Av,
                 'object_dist': object_dist,
                 'object_sfr_mult_factorA': 1e9,
                 'file_imf': get_imf(gal.target),
                 'binary_kind': 1,
                 'binary_frac': 0.35}

    gal_inp.add_params(galinp_kw)

    galaxy_input = os.path.join(snap_src, 'input', 'input_%s.dat' % gal.target)
    gal_inp.write_params(galaxy_input, rsp.TrilegalUtils.galaxy_input_fmt())
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
                                      count_offset=count_offset, photsys=gal.photsys)
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
                                              count_offset=count_offset, photsys=gal.photsys)

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
    
    sgal = rsp.Galaxies.simgalaxy(sim_file, filter1=gal.filter1, filter2=gal.filter2,
                                  photsys=gal.photsys)
    sgal.load_ast_corrections()

    sgal = make_normalized_simulation(gal, model, 'gal', 'gal', maglims=maglims,
                                      band=band, offsets=offsets,
                                      run_trilegal=False, leo_norm=leo_norm,
                                      trilegal_output=sim_file)

    return gal, sgal


def call_mc_norm(targets=None, models=None, gi10=False, **mc_norm_kw):
    '''
    calls mc_norm for many targets, over many models. 
    '''

    targets = load_targets(targets)

    if models is None:
        models = ['cmd_input_CAF09_S_JAN13.dat',
                  'cmd_input_gi10_rev_old_tracks.dat',
                  'cmd_input_gi10_rev.dat']
    if type(models) == str:
        models = [models]
        
    mc_norm_kw = dict({'band': 'ir', 'offsets': (1.5, 0.), 'leo_norm': False,
                       'maglims': 'trgb'}.items() + mc_norm_kw.items())
    
    if gi10 is True:
        targets = gi10_overlap()
        models = ['cmd_input_gi10_bow.dat']
        mc_norm_kw = {'band': 'opt', 'offsets': (2., 0.), 'leo_norm': True,
                      'maglims': 'trgb', 'leo_method': True}

    for target in targets:
        for model in models:
            mc_norm(target, model, **mc_norm_kw)
            #sgal_rgb_agb(target, model, **mc_norm_kw)

def mc_norm(target, model, band=None, input_file=None, maglims=None,
            offsets=None, leo_norm=False, leo_method=False, **kwargs):
    '''
    One random draw to make the LF plots could be any number of TP-AGB stars
    within the Poisson noise, or perhaps even wider range. This code does the
    normalization 10,001 times to see what the mean Nagb/Nrgb ratio is and its
    stdev.
    '''    
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
        pct_diff = (np.mean(nmodel) - nrgb_nagb_data) / nrgb_nagb_data

    print '%s %.2f %.2f %.4f %.2f %s' % (target, nrgb_nagb_data, np.mean(nmodel),
                                         np.std(nmodel), pct_diff, model)


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


def targets_z002():
    return ['DDO82', 'IC2574-SGS', 'UGC-4305-1', 'UGC-4305-2', 'NGC4163', 'UGC8508']

def targets_paper1():
    return np.concatenate([gi10_overlap(), targets_z002()])
    #return ['UGC-4305-1', 'UGC-4305-2', 'NGC4163', 'UGC8508']

def targets_leo_norm_ok():
    '''
    targets that don't have MS to cut out, might as well use the whole
    RGB
    '''
    return np.concatenate([gi10_overlap(), ['DDO82', 'NGC4163', 'UGC8508']])
    
def load_ast_file(gal, model):
    model_short = model.replace('cmd_input_','').lower()
    search_str = '*'.join(('ast_%s' % gal.filter1, gal.filter2, gal.target,
                           model_short))
    try:
        sim_file, = rsp.fileIO.get_files(os.path.join(snap_src, 'output'),
                                         search_str)
    except ValueError:
        raise ValueError, 'more than one or zero value(s) found when searching %s' % search_str
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


def load_targets(targets):
    if targets == 'all':
        targets = all_targets()
    elif targets == 'gi10':
        targets = gi10_overlap()
    elif targets == 'z002':
        targets = targets_z002()
    elif targets == 'paper1':
        targets = targets_paper1()

    if type(targets) == str:
        targets = [targets]
    return targets


def nir_cmd_plot_limits(xlim=None, ylim=None, xlim2=None):
    # these are chosen for IR plots.
    if xlim is None:
        xlim = (-0.5, 1.5)
    if ylim is None:
        ylim = (25.3, 18)
    if xlim2 is None:
        xlim2 = (0.8, 4e4)
    return xlim, ylim, xlim2 


def ir_rgb_agb_ratio(renormalize=False, filt1=None, filt2=None, band=None,
                     targets=None, run_trilegal=False, leo_ast=True, 
                     offsets=(1.5, 0.), models=None, maglims=None,
                     leo_method=False, leo_norm=False, make_plot=True,
                     xlim=None, ylim=None, xlim2=None, add_boxes=True,
                     color_hist=False,  plt_tpagb=False, **kwargs):
    smgs = []
    targets = load_targets(targets)

    if type(models) == str:
        models = [models]

    norm_sim_kw = {'offsets': offsets, 'band': band, 'leo_ast': leo_ast,
                   'run_trilegal': run_trilegal, 'leo_method': leo_method,
                   'leo_norm': leo_norm}

    if make_plot is True:
        xlim, ylim, xlim2 = nir_cmd_plot_limits(xlim=xlim,
                                                ylim=ylim,
                                                xlim2=xlim2)

        plot_LF_kw = {'ylim': ylim, 'xlim': xlim, 'xlim2': xlim2,
                      'color_hist': color_hist, 'title': True}

    result_dict = {}
    for target in targets:
        logger.info('working on %s' % target)
        norm_sim_kw['object_mass'] = load_sim_masses(target)
        gal = load_galaxy(target, band=band)
        leo_now = norm_sim_kw['leo_method']
        if target in targets_leo_norm_ok():
            norm_sim_kw['leo_method'] = True
        else:
            norm_sim_kw['leo_method'] = False
        if leo_now != norm_sim_kw['leo_method']:
            logger.info('over-ridcing input, leo_method is now %s' % leo_method)

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

            pct_diff = (nrgb_nagb_sim - nrgb_nagb_data) / nrgb_nagb_data

            if make_plot is True:
                plot_LF_kw['model_title'] = translate_model_name(model_title)
                smg.make_LF(gal.filter1, gal.filter2, plot_LF_kw=plot_LF_kw,
                            comp50=True, add_boxes=add_boxes,
                            color_hist=color_hist, plot_tpagb=plot_tpagb)
            
            sub_dict = {'data': nrgb_nagb_data, model_title: nrgb_nagb_sim,
                        'pct_diff': pct_diff}

            if result_dict.has_key(target):
                result_dict[target].update(sub_dict)
            else:
                result_dict[target] = sub_dict
            smgs.append(smg)
    #print '   ', ' '.join(result_dict[result_dict.keys()[0]].keys())
    #for key, rd in result_dict.items():
    #    print key, ' '.join(['%.3f' % i for i in  rd.values()])
    pprint(result_dict)
    return smgs

def match_tests(target, model, band, verts=False, inverse_verts=False, inputs={},
                extra='', match_phot=None):
    '''
    run the match fitter with the simulated galaxy as the bg,
    no time bins, and data as data. this will just do the chi2 test.
    
    verts will load the agb verts and only compare that area of a hess
    diagram.
    
    inverse_verts will ignore the agb verts section of the cmd.
    '''
    
    gal, sgal = load_normalized_simulation(target, model, band=band, 
                                           maglims=inputs['maglims'],
                                           offsets=inputs['offsets'],
                                           leo_norm=inputs['leo_norm'])
    models_loc = '/Users/phil/research/TP-AGBcalib/SNAP/models/match/'
    fmter = '_'.join(['match', extra, gal.target, gal.filter1, gal.filter2, model]).replace('.dat', '').lower()

    match_bg = os.path.join(models_loc, 'bg', fmter + '.bg')

    color = sgal.ast_color[sgal.norm_inds]
    mag2 = sgal.ast_mag2[sgal.norm_inds]

    if match_phot is None:
        match_phot, = rsp.fileIO.get_files(os.path.join(models_loc, 'phot'),
                                           '*%s*.gst.match' % gal.target)

    if verts is True or inverse_verts is True:
        agb_verts = load_agb_verts(gal)

        # load photometry verts, will need to save new file
        m1, m2 = np.loadtxt(match_phot, unpack=True)
        phot_points = np.column_stack((m1-m2, m2))
        phot_mask = nxutils.points_inside_poly(phot_points, agb_verts)
        
        # sim cmd points
        points = np.column_stack((color, mag2))
        mask = nxutils.points_inside_poly(points, agb_verts)

        if verts is True:
            vinds, = np.nonzero(mask)
            vpinds, = np.nonzero(phot_mask)
            match_phot_ext = '_verts.dat'
        if inverse_verts is True:
            vinds, = np.nonzero(mask==0)
            vpinds, = np.nonzero(phot_mask==0)
            match_phot_ext = '_inverts.dat'

        photpts = len(vpinds)
        simpts = len(vinds)

        # slice sim cmd
        color = color[vinds]
        mag2 = mag2[vinds]
        # save new match phot
        match_phot = rsp.fileIO.replace_ext(match_phot, match_phot_ext)
        np.savetxt(match_phot, np.array([m1[vpinds], m2[vpinds]]).T)
    
    photpts = len(gal.mag2)
    simpts = len(mag2)

    rsp.match_utils.write_match_bg(color, mag2, match_bg)

    pm_file = os.path.join(models_loc, 'pars', fmter + '.par')

    match_fake = get_fake_files(gal.target, band=band)
    match_out = os.path.join(models_loc, 'output', fmter + '.fit')
    msg = os.path.join(models_loc, 'msgs', fmter + '.msg')
    figname = os.path.join(models_loc, 'plots', fmter + '.png')

    match_kwargs = {'dmod': sgal.data.get_col('m-M0')[0],
                    'Av': sgal.data.get_col('Av')[0],
                    'filter1': sgal.filter1.replace('F', 'IR'),
                    'filter2': sgal.filter2.replace('F', 'IR'),
                    'pmfile': pm_file,
                    'color': color,
                    'mag2': mag2,
                    'match_bg': match_bg#,
                    #'bright1': agb_verts[:, 1].min() + agb_verts[:, 0].min(),
                    #'faint1': agb_verts[:, 1].max() + agb_verts[:, 0].max(),
                    #'bright2': agb_verts[:, 1].min(),
                    #'faint2': agb_verts[:, 1].max(),
                    #'colmin': agb_verts[:, 0].min(),
                    #'colmax': agb_verts[:, 0].max()
                    }

    npts = np.min([photpts, simpts])
    if verts is True:
        dmag = np.diff(np.linspace(m2.min(), m2.max(), int(np.sqrt(npts))))[0]
        dcol = np.diff(np.linspace(np.min((m1-m2)), np.max((m1-m2)), int(np.sqrt(npts))))[0]
        dmag = np.min([0.15, dmag])
        dcol = np.min([0.1, dcol])
        match_kwargs['dmag'] = dmag
        match_kwargs['dcol'] = dcol
        print match_kwargs['dmag'], match_kwargs['dcol']
    
    model_name = translate_model_name(model)

    chi2, fit = rsp.match_utils.match_light(gal, pm_file, match_phot,
                                            match_fake, match_out, msg,
                                            match_kwargs=match_kwargs,
                                            make_plot=True, model_name=model_name,
                                            figname=figname, loud=True)
    return chi2, fit

def translate_model_name(model):
    if 'apr13vw93' in model.lower():
        model_name = '$\dot{M}_{VW93}$'
    elif 'apr13' in model.lower():
        model_name = '$\dot{M}_{BS95}$'
    elif 'gi10' in model.lower() or 'bow' in model.lower():
        model_name = '$\dot{M}_{G10}$'
    elif 'mar13' in model.lower():
        model_name = '$\dot{M}_{M13}$'
    return model_name

def fuck_you_match(gal, sgal, verts=False, inverse_verts=False, make_plot=False,
                   figname=None, dmag=0.1, dcol=0.05, modelB=None,
                   rgb_only=False):

    if rgb_only is True:
        sgal.all_stages('RGB')
        inds = np.intersect1d(sgal.irgb, sgal.norm_inds)
    else:
        inds = sgal.norm_inds
    sim_color = sgal.ast_color[inds]
    sim_mag2 = sgal.ast_mag2[inds]

    if modelB is None:
        data_color = gal.color
        data_mag2 = gal.mag2
    else:
        if rgb_only is True:
            gal.all_stages('RGB')
            ginds = np.intersect1d(gal.irgb, gal.norm_inds)
        else:
            ginds = gal.norm_inds

        data_color = gal.ast_color[ginds]
        data_mag2 = gal.ast_mag2[ginds]

    if verts is True or inverse_verts is True:
        agb_verts = load_agb_verts(gal)

        # load photometry verts, will need to save new file
        phot_points = np.column_stack((data_color, data_mag2))
        phot_mask = nxutils.points_inside_poly(phot_points, agb_verts)

        # sim cmd points
        points = np.column_stack((sim_color, sim_mag2))
        mask = nxutils.points_inside_poly(points, agb_verts)

        if verts is True:
            vinds, = np.nonzero(mask)
            vpinds, = np.nonzero(phot_mask)
            match_phot_ext = '_verts.dat'
        if inverse_verts is True:
            vinds, = np.nonzero(mask==0)
            vpinds, = np.nonzero(phot_mask==0)
            match_phot_ext = '_inverts.dat'

        photpts = len(vpinds)
        simpts = len(vinds)

        # slice sim cmd
        sim_color = sim_color[vinds]
        sim_mag2 = sim_mag2[vinds]
        data_color = data_color[vpinds]
        data_mag2 = data_mag2[vpinds]
        
    # spacing
    photpts = len(data_mag2)
    simpts = len(sim_mag2)
    npts = np.min([photpts, simpts])
    dmag = np.diff(np.linspace(np.min(data_mag2), np.max(data_mag2),
                               int(np.sqrt(npts))))[0]
    dcol = np.diff(np.linspace(np.min(data_color), np.max(data_color),
                               int(np.sqrt(npts))))[0]
    dmag = np.min([0.15, dmag])
    dcol = np.min([0.1, dcol])

    data_hess = rsp.astronomy_utils.hess(data_color, data_mag2, dmag,
                                         **{'cbinsize': dcol})

    sim_hess = rsp.astronomy_utils.hess(sim_color, sim_mag2, dmag,
                                         **{'cbin': data_hess[0],
                                            'mbin': data_hess[1]})

    chi2, pct_dif, sig = rsp.Galaxies.stellar_prob(data_hess, sim_hess)


    if make_plot is True:
        if figname is None:
            figname = 'pgpro_%s_%s_%s_%s%s.png' % (gal.target, model, gal.filter1, gal.filter2, extra)        

        model_name = translate_model_name(sgal.model)

        ZS = [data_hess[2], sim_hess[2], pct_dif, sig]
        xlim, ylim, _ = nir_cmd_plot_limits()

        extent = [data_hess[0][0], data_hess[0][-1], data_hess[1][-1], data_hess[1][0]]
        labels = ['${\\rm %s}$' % gal.target, model_name, '$\%\ {\\rm Difference}$', '$\chi^2=%.2f$' % chi2]
        grid = rsp.match_graphics.match_plot(ZS, extent, xlim, ylim,
                                             labels=labels,
                                             **{'xlabel': '$%s-%s$' % (gal.filter1, gal.filter2),
                                                'ylabel': '$%s$' % gal.filter2})
        plt.savefig(figname, dpi=300)
        plt.close()

    return chi2


def get_nrgb_from_optical(targets, trgb_offset=1.5):
    Nrgbs = []
    targets = load_targets(targets)
    for target in targets:
        gal = load_galaxy(target, band='opt')

        if target in gi10_overlap():
            cmin, cmax = -100, 100
        elif gal.filter1 == 'F475W':
            cmin, cmax = 1.5, 2.6
        elif gal.filter1 == 'F606W':
            cmin, cmax = 0.5, 1.75
        elif gal.filter1 == 'F555W':
            cmin, cmax = 0.8, 2.0
        else:
            print 'what the fuck are you trying to do here?'
    
        mag_faint, mag_bright = gal.Trgb + trgb_offset, gal.Trgb
    
        big_box = np.array([[cmin, mag_faint],
                            [cmin, mag_bright],
                            [cmax, mag_bright],
                            [cmax, mag_faint],
                            [cmin, mag_faint]])
        points = np.column_stack((gal.Color, gal.Mag2))
        inside, = np.nonzero(nxutils.points_inside_poly(points, big_box))
        
        if target in gi10_overlap():
            inds = inside
        else:
            inds = inside

        #    # the 1/2 std around color mean
        #    cmean = np.mean(gal.Color[inside])
        #    cstd = np.std(gal.Color[inside])/2.
        #    col_min = cmean - cstd
        #    col_max = cmean + cstd
        #    rgb_box = np.array([[col_min, mag_faint],
        #                        [col_min, mag_bright],
        #                        [col_max, mag_bright],
        #                        [col_max, mag_faint],
        #                        [col_min, mag_faint]])
        #
        #    inds, = np.nonzero(nxutils.points_inside_poly(points, rgb_box))    
        #
        Nrgb = len(inds)
        print gal.target, gal.filter1, gal.filter2, Nrgb
        Nrgbs.append(Nrgb)
        fig, ax = gal.plot_cmd(gal.Color, gal.Mag2, scatter_off=True)
        ax.plot(gal.Color[inds], gal.Mag2[inds], '.')
        ax.set_title(gal.target)
        if not target in gi10_overlap():
            from matplotlib.pyplot import ginput
            pts = ginput(0, timeout=-1)
            xs, ys = np.transpose(pts)
            inds, = np.nonzero(nxutils.points_inside_poly(points, pts))
            print len(inds)
            print pts
    return Nrgbs

def run_fuck_you_match(targets=None, models=None, band=None, inputs={},
                       modelB=None):
    '''
    does the chi2 test for full field, agb region, and not-agb region.
    to compare model to model, use modelB.
    '''

    if targets is None:
        targets = inputs['targets']
    if band is None:
        band = inputs['band']
    if models is None:
        models = inputs['models']

    targets = load_targets(targets)

    if type(models) == str:
        models = [models]

    figname_fmt = 'pgpro_%s_%s_%s.png'

    for target in targets:
        print target
        for model in models:
            print model
            gal, sgal = load_normalized_simulation(target, model, band=band, 
                                                   maglims=inputs['maglims'],
                                                   offsets=inputs['offsets'],
                                                   leo_norm=inputs['leo_norm'])

            if modelB is not None:
                _, sgal2 = load_normalized_simulation(target, modelB, band=band, 
                                                      maglims=inputs['maglims'],
                                                      offsets=inputs['offsets'],
                                                      leo_norm=inputs['leo_norm'])

                # overwrite gal with modelB.
                [gal.__setattr__(k, v) for k,v in sgal2.__dict__.items()]
                #gal.target = modelB.replace('.dat', '').split('_')[-1]
            # entire cmd fit
            extra = 'full'
            print extra
            figname = figname_fmt % (target, model.replace('.dat',''), extra)
            full_chi2 = fuck_you_match(gal, sgal, figname=figname, modelB=modelB,
                                       make_plot=True)
            
            if band == 'opt':
                print '%s %s %.2f' % (target, model.replace('.dat', '').split('_')[-1], full_chi2)
                return
            # include agb only
            extra = 'agb'
            print extra
            figname = figname_fmt % (target, model.replace('.dat',''), extra)
            agb_chi2 = fuck_you_match(gal, sgal, figname=figname, verts=True,
                                      make_plot=True, modelB=modelB)
            # exclude agb
            extra = 'not_agb'
            print extra
            figname = figname_fmt % (target, model.replace('.dat',''), extra)
            part_chi2 = fuck_you_match(gal, sgal, figname=figname, modelB=modelB,
                                       inverse_verts=True, make_plot=True)

            print '%s %s %.2f %.2f %.2f' % (target, model.replace('.dat', '').split('_')[-1], full_chi2, agb_chi2, part_chi2)


def run_match_tests(targets=None, models=None, band=None, inputs={},
                    modelB=None, ex2=''):
    '''
    calls match_tests in a batch mode.
    '''

    if targets is None:
        targets = inputs['targets']
    if band is None:
        band = inputs['band']
    if models is None:
        models = inputs['models']

    targets = load_targets(targets)

    if type(models) == str:
        models = [models]

    if modelB is not None:
        bg_loc = '/Users/phil/research/TP-AGBcalib/SNAP/models/match/bg/'
        match_phot_fmt = 'match_%s_%s*%s.bg'

    for target in targets:
        print target
        for model in models:
            print model
            # entire cmd fit
            extra = 'full'
            match_phot, = rsp.fileIO.get_files(bg_loc, str.lower(match_phot_fmt % (extra, target, modelB.replace('.dat',''))))
            full_chi2, full_fit = match_tests(target, model, band, inputs=inputs,
                                              extra=extra+ex2, match_phot=match_phot)
            # include agb only
            extra = 'agb'
            match_phot, = rsp.fileIO.get_files(bg_loc, str.lower(match_phot_fmt % (extra, target, modelB.replace('.dat',''))))
            agb_chi2, agb_fit = match_tests(target, model, band, verts=True,
                                            inputs=inputs, extra=extra+ex2,
                                            match_phot=match_phot)
            # exclude agb
            extra = 'not_agb'
            match_phot, = rsp.fileIO.get_files(bg_loc, str.lower(match_phot_fmt % (extra, target, modelB.replace('.dat',''))))
            part_chi2, part_fit = match_tests(target, model, band,
                                              inverse_verts=True, inputs=inputs,
                                              extra=extra+ex2,
                                              match_phot=match_phot)

            print target, model, full_chi2, agb_chi2, part_chi2
            #print target, model, agb_chi2

def check_stats():
    filename = '/Users/phil/Desktop/tpagb_stats.dat'
    lines = open(filename, 'r').readlines()
    target, data, model, model_std, pct_diff, model_name = zip(*[l.strip().split() for l in lines])
    data = np.array(data, dtype=float)
    model = np.array(model, dtype=float)
    model_std = np.array(model_std, dtype=float)
    pct_diff = np.array(pct_diff, dtype=float)
    target = np.array(target, dtype=str)
    model_name = np.array(model_name, dtype=str)

    isort = np.argsort(data)
    fig, ax = plt.subplots()
    for i in range(4):
        ax.plot(pct_diff[isort[i::4]], data[isort[i::4]], label=model_name[i])
    ax.legend()

    filename = '/Users/phil/Desktop/chi2tests.dat'
    lines = open(filename, 'r').readlines()
    target, model, full, agb, not_agb = zip(*[l.strip().split() for l in lines if not l.startswith('#')])
    model = np.array(model, dtype=str)
    model[model=='APR13'] = '$BS95$'
    model[model=='bow'] = '$G10$'
    model[model=='MAR13'] = '$M13$'
    model[model=='APR13VW93'] = '$VW93$'
    
    full = np.array(full, dtype=float)
    target = np.array(['$%s$' % t for t in target], dtype=str)
    agb = np.array(agb, dtype=float)
    not_agb = np.array(not_agb, dtype=float)

    # a plot of the chi full minus chi no agb.
    cols = ['darkred', 'navy', 'black', 'purple']
    fig, ax = plt.subplots()
    for i in range(4):
        ax.plot(full[i::4] - not_agb[i::4], 'o', ms=15, alpha=0.4, lw=2, color=cols[i], label=model[i])
        #ax.plot(not_agb[i::4], 'o', lw=2, , color=cols[i])
    ax.legend(loc=0, numpoints=1)
    plt.xticks(np.arange(len(full[i::4])), target[i::4], rotation=45)
    ax.set_xlim(-0.5, 9.5)
    ax.set_ylabel('$\chi^2_{full}-\chi^2_{no-agb}$', fontsize=20)
    plt.subplots_adjust(bottom=0.2)
    plt.savefig('chi2_agb_improves_fit.png', dpi=300)

    # a plot of just chi2 with full and no agb.
    cols = ['darkred', 'navy', 'black', 'purple']
    fig, ax = plt.subplots()
    for i in range(4):
        ax.plot(full[i::4], 'o', ms=15, alpha=0.4, lw=2, color=cols[i], label=model[i])
        ax.plot(not_agb[i::4], 'x', ms=15, lw=2, color=cols[i])
    ax.legend(loc=0, numpoints=1)
    plt.xticks(np.arange(len(full[i::4])), target[i::4], rotation=45)
    ax.set_xlim(-0.5, 9.5)
    ax.set_ylabel('$\chi^2$', fontsize=20)
    plt.subplots_adjust(bottom=0.2)
    plt.savefig('chi2_full_no_agb.png', dpi=300)


    model = model.reshape(9,4)        
    not_agb = not_agb.reshape(9,4)        
    full = full.reshape(9,4)        
    cols = ['darkred', 'navy', 'purple', 'purple']
    fig, ax = plt.subplots()
    for i in range(len(model)):
        # G10 - x: 
        fg10 = full[i][2]
        ng10 = not_agb[i][2]
        if np.isnan(ng10):
            print 'nan!!'
            ng10 = 1e10
        full_diff = np.array([fg10 - full[i][j] for j in [0,1,3]])
        nagb_diff = np.array([ng10 - not_agb[i][j] for j in [0,1,3]])
        #nagb_diff = not_agb[i][2] - not_agb[i]
        [ax.plot(i, full_diff[j], 'o', ms=10, alpha=0.4, color=cols[j]) for j in range(len(full_diff))]
        [ax.plot(i, nagb_diff[j], 'x', ms=10, color=cols[j]) for j in range(len(nagb_diff))]
    [ax.plot(-999, -999, 'o', ms=10, alpha=0.4, color=cols[j], label=model[i][j]) for j in [0,1,3]]
    ax.set_xlim(-0.5, 9.5)
    ax.hlines(0, *ax.get_xlim(), linestyle='--')
    #ax.set_ylim(-.2, 0.7)
    ax.legend(loc=0, numpoints=1)
    plt.xticks(np.arange(len(target[::4])), target[::4], rotation=45)
    plt.subplots_adjust(bottom=0.2)
    ax.set_ylabel('$\chi^2_{G10}-\chi^2_i$', fontsize=20)
    plt.savefig('chi2_gi10_comp.png', dpi=300)
    # left: 12
    # bottom: 22
    # right: 96
    # top: 96 
    # w,hspace: .2 


def compare_LFs(smgs):

    bmap = brewer2mpl.get_map('Paired', 'Qualitative', 9)
    cols = bmap.mpl_colors

    models, imod = np.unique([smg.sgal.model for smg in smgs], return_inverse=True)

    # j will go in order of galaxy.
    for i in range(len(models)):
        resids = []
        pct_diffs = []
        fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, figsize=(10,12))
        model_name = translate_model_name(models[i])
        for j, smg in enumerate(np.array(smgs)[np.where(imod==i)]):
            data_hist = smg.gal_hist
            model_hist = smg.sgal_hist
            bins = rsp.astronomy_utils.mag2Mag(smg.bins[1:], 'F160W',
                                               'wfc3snap',
                                               **{'target': smg.gal.target,
                                                  'filter1': smg.gal.filter1})
            #bins = smg.bins[1:]
            comp_inds = np.nonzero(bins < smg.gal.trgb)
            pct_diff = (model_hist - data_hist) / data_hist
            resid = model_hist - data_hist
            resids.append(np.mean(resid[comp_inds]))
            pct_diffs.append(np.mean(pct_diff[comp_inds]))
            # 3 panel
            # top: # vs magnitude
            # middle: residual model-data  vs mag
            # bottom: %diff vs mag
            label='$%s$' % smg.gal.target.replace('-','\!-\!')
            plt_kw = {'ls': 'steps', 'color': cols[j], 'lw': 2}
            #[ax.vlines(smg.gal.Trgb, *ax.get_ylim(), color=cols[j]) for ax in [ax1, ax2, ax3]]
            ax1.semilogy(bins[comp_inds], model_hist[comp_inds], label=label, **plt_kw)
            ax2.plot(bins[comp_inds], pct_diff[comp_inds], **plt_kw)
            ax3.plot(bins[comp_inds], resid[comp_inds], **plt_kw)
            plt_kw['ls'] += '--'
            plt_kw['lw'] = 1
            ax1.semilogy(bins[comp_inds], data_hist[comp_inds], **plt_kw)
            print smg.gal.target, np.mean(resid[comp_inds]), np.std(resid[comp_inds]), np.mean(pct_diff[comp_inds]), np.std(pct_diff[comp_inds])
        print model_name, np.mean(resids), np.std(resids), np.mean(pct_diffs), np.std(pct_diffs)
        ax1.text(0.5, 0.85, '$%s$' % model_name, fontsize=20, transform=ax1.transAxes)
        ax2.text(0.1, 0.85, '${\\rm \%\ Difference}$', fontsize=20,
                 transform=ax2.transAxes)        
        ax3.text(0.1, 0.85, '${\\rm Residual}$', fontsize=20,
                 transform=ax3.transAxes)

        ax1.legend(loc=0, numpoints=1, frameon=False)
        ax1.set_ylabel('$\#$', fontsize=20)
        ax2.set_ylabel('$(m-n)/n$', fontsize=20)
        ax3.set_ylabel('$m-n$', fontsize=20)
        [ax.set_xlabel('$F160W$', fontsize=20) for ax in [ax1, ax2, ax3]]
        plt.tick_params(labelsize=16)
        ax1.set_ylim(0, 400)
        ax2.set_ylim(-1.5, 30)
        ax3.set_ylim(-75, 275)
        [ax.set_xlim(-10, -5.5) for ax in [ax1, ax2, ax3]]
        plt.subplots_adjust(hspace=.25, top=.95)
        plt.savefig('comp_lfs_%s.png' % models[i].replace('.dat','').split('_')[-1], dpi=300)
    #fig, (axs) = plt.subplots(figsize=(10,10), nrows=2, ncols=2)
    #axs = axs.ravel()
    do_a_plot_that_wastes_time = True
    if do_a_plot_that_wastes_time is True:
        for i in range(len(models)):
            model_name = translate_model_name(models[i])
            fig, ax  = plt.subplots(figsize=(10, 10))

            divider = make_axes_locatable(ax)
            axt = divider.append_axes("top", 1.5, pad=0.1, sharex=ax)
            #axl = divider.append_axes("left", 1.5, pad=0.1, sharey=ax)
            axr = divider.append_axes("right", 1.5, pad=0.1, sharey=ax)
            # this is for an awesome four panel that I spent to much time
            # on and feel stupid deleting.
            #divider.new_vertical(1.5, pad=0.1, pack_start=True)
            #plt.draw()
            #mid = ax.get_position().get_points()
            #top = axt.get_position().get_points()
            #bot = (mid[0][0], 0.1, top[1][0]-top[0][0],  top[1][1]-top[0][1])
            #axb = plt.axes(bot, sharex=ax)

            ax.text(1.12, 1.12, model_name, fontsize=20, transform=ax.transAxes)
            for j, smg in enumerate(np.array(smgs)[np.where(imod==i)]):
                plt_kw = {'histtype': 'step', 'color': cols[j], 'lw': 2}
                smg.sgal.all_stages('TPAGB')
                
                #smg.sgal.convert_mag(target=smg.gal.target)
                points = np.column_stack((smg.sgal.ast_color, smg.sgal.ast_mag2))
                inregion, = np.nonzero(nxutils.points_inside_poly(points, smg.agb_verts))

                tpagb_inds = list(set(smg.sgal.norm_inds) & set(smg.sgal.itpagb) & set(inregion))
                poop = np.intersect1d(smg.sgal.norm_inds, inregion) 
                not_tpagb = np.array([k for k in poop if k not in smg.sgal.itpagb])
                label='$%s$' % smg.gal.target.replace('-','\!-\!')

                mag2 = rsp.astronomy_utils.mag2Mag(smg.sgal.ast_mag2, 'F160W',
                                               'wfc3snap',
                                               **{'target': smg.gal.target,
                                                  'filter1': smg.gal.filter1})
                x = mag2[tpagb_inds]
                xx = mag2[not_tpagb]
                ydata = smg.sgal.data.get_col('logAge')
                #ydata = smg.sgal.data.get_col('m_ini')
                y = ydata[tpagb_inds]
                yy = ydata[not_tpagb]
    

                #ax.plot(xx, yy, 'o', ms=3, alpha=.5, color=cols[j], mec=cols[j])
                ax.plot(x, y, 'o', ms=5, alpha=.5, color=cols[j], mec=cols[j], label=label)
            
                ax.set_xlim(-10., -5.5)
                #ax.set_ylim(0.5, 6.)
                ax.set_ylim(10.5, 7.5)
                # make some labels invisible
                plt.setp(axt.get_xticklabels() + axr.get_yticklabels(),
                         visible=False)

                bins = int(np.round(np.sqrt(len(x))))
                axt.hist(x, bins=bins, **plt_kw)
                axr.hist(y, bins=bins, log=True, orientation='horizontal', **plt_kw)
                #axb.hist(xx, bins=bins, log=True, **plt_kw)
                #axl.hist(yy, bins=bins, orientation='horizontal', **plt_kw)
            #axl.set_xlim(axl.get_xlim()[::-1])
            #axl.xaxis.set_major_locator(MaxNLocator(4))
            #axr.xaxis.set_major_locator(MaxNLocator(4))
            #axb.xaxis.set_major_locator(MaxNLocator(4))
            #axb.xaxis.set_major_locator(MaxNLocator(6))
            
            ax.set_xlabel('$F160W$', fontsize=20)
            #ax.set_ylabel('$M (M_\odot)$', fontsize=20)
            ax.set_ylabel('${\log \\rm Age\ (yr)}$', fontsize=20)
            #axt.text(.1,.85, '$TP-AGB$', transform=axt.transAxes)
            #axr.text(.1,.95, '$TP-AGB$', transform=axr.transAxes)
            #axb.text(.1,.85, '${\rm All\ Stars}$', transform=axb.transAxes)
            for tl in axt.get_xticklabels():
                tl.set_visible(False)

            for tl in axr.get_yticklabels():
                tl.set_visible(False)
            #axHisty.set_xticks([0, 50, 100])
            #[ax.set_ylim(-2, -10) for ax in axs]
            #[ax.set_xlim(0, 6) for ax in axs]
            ax.legend(loc=0, numpoints=1, frameon=False)
            plt.draw()
            plt.savefig('age_mag_%s.png' % models[i].replace('.dat','').split('_')[-1], dpi=300)
        #[ax.vlines(smg.gal.Trgb, *ax.get_ylim(), color=cols[j]) for ax in [ax1, ax2, ax3]]

def match_metals():
    '''
    print the sfh weighted metallicity from the match input file. 
    '''
    sfh_loc = '/Users/phil/research/TP-AGBcalib/SNAP/data/sfh/'
    sfh_files = rsp.fileIO.get_files(sfh_loc, '*.dat')
    for sfh_file in sfh_files:
        a, s, z = np.loadtxt(sfh_file, unpack=True)
        o, = np.nonzero(np.log10(a*1e9) < 8.5)
        print '%.4f %.4f %.4f %s ' % (np.min(z), np.sum(z*s)/np.sum(s), np.max(z), sfh_file.split('_')[2])
        #print '%.4f %.4f %.4f %s ' % (np.min(z[o]), np.sum(z[o]*s[o])/np.sum(s[o]), np.max(z[o]), sfh_file.split('_')[2])

def agb_logl_age():
    for i in range(len(models)):
        model_name = models[i].replace('.dat', '').split('_')[-1]
        agb_track_loc = '/Users/phil/research/TP-AGBcalib/AGBTracks/CAF09/S_%s/S12_Z0.002_Y0.252/' % model_name
        if not os.path.isdir(agb_track_loc) is True:
            print model_name, 'no agb tracks found'
        model_name = translate_model_name(models[i])
        agb_track_names = [os.path.join(agb_track_loc, a)
                           for a in os.listdir(agb_track_loc)
                           if a.startswith('agb_')]
        tracks = [fileIO.get_numeric_data(agb_track)
                  for agb_track in agb_track_names]

        logls = np.array([t.get_col('L_star') for t in tracks])
        logts = np.array([t.get_col('T_star') for t in tracks])
        ages =  np.array([t.get_col('ageyr') for t in tracks])
        [ax.plot(ages[i], logls[i]) for i in range(len(ages))]
    
def agb_lifetimes(models):
    fig, ax = plt.subplots()
    for i in range(len(models)):
        model_name = models[i].replace('.dat', '').split('_')[-1]
        agb_track_loc = '/Users/phil/research/TP-AGBcalib/AGBTracks/CAF09/S_%s/S12_Z0.002_Y0.252/' % model_name
        if not os.path.isdir(agb_track_loc) is True:
            print model_name, 'no agb tracks found'
        model_name = translate_model_name(models[i])
        agb_track_names = [os.path.join(agb_track_loc, a)
                           for a in os.listdir(agb_track_loc)
                           if a.startswith('agb_')]
        tracks = [fileIO.get_numeric_data(agb_track)
                  for agb_track in agb_track_names]
        masses = np.array([t.mass for t in tracks])
        statuss = t.data_array['status']
        taus = np.array([np.sum(t.data_array['dt']) for t in tracks])
        ax.plot(masses, taus/1e6, lw=2, label=model_name)

    ax.set_xlabel('${\\rm Initial\ Mass\ (M_\odot)}$', fontsize=20)
    ax.set_ylabel('${\\rm Lifetime\ (Myr)}$', fontsize=20)
def main(inputfile):
    #ch = logging.StreamHandler()
    #ch.setLevel(logging.DEBUG)
    #logger.addHandler(ch)

    #import cProfile
    #command = 'ir_rgb_agb_ratio()'
    #cProfile.runctx(command, globals(), locals())
    inputs = rsp.fileIO.load_input(inputfile)
    if inputs['debug'] is True:
        import pdb
        pdb.set_trace()
        del inputs['debug']

    global snap_src
    snap_src = inputs['snap_src']
    del inputs['snap_src']

    mc_norm = inputs['mc_norm']
    del inputs['mc_norm']
    
    if inputs['make_lfs'] is True:
        ir_rgb_agb_ratio(**inputs)
    
    if mc_norm is True:
        call_mc_norm(**inputs)

    if inputs['stats_tests'] is True:
        run_fuck_you_match(inputs=inputs)

if __name__ == "__main__":
    main(sys.argv[1])
else:
    global snap_src
    snap_src = '/Users/phil/research/TP-AGBcalib/SNAP'