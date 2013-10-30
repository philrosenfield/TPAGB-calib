import logging
logging.basicConfig(filename='galaxy_tests.log',level=logging.DEBUG)
logger = logging.getLogger()
logger.info('start of run')
import ResolvedStellarPops as rsp
import os
import sys
import difflib
import numpy as np

import matplotlib.pyplot as plt
import matplotlib.nxutils as nxutils
from mpl_toolkits.axes_grid1 import make_axes_locatable
from matplotlib.colors import LogNorm
from pprint import pprint
import brewer2mpl
from TPAGBparams import research_path, table_src
import multiprocessing
angst_data = rsp.angst_tables.AngstTables()


class VarySFHs(rsp.match_utils.StarFormationHistories):
    def __init__(self, galaxy_input, match_sfh_file, outfile_loc='default'):
        super(VarySFHs, self).__init__(match_sfh_file)
        self.galaxy_input = galaxy_input
        if outfile_loc == 'default':
            self.outfile_loc = '/home/rosenfield/research/TP-AGBcalib/SNAP/data/sfh_parsec/templates/ddo71/mc/inputs/'
        rsp.fileIO.ensure_dir(self.outfile_loc)

    def prepare_trilegal_sfr(self, make_many_kw=None):
        make_many_kw = make_many_kw or {}
        self.sfr_files = self.make_many_trilegal_sfhs(**make_many_kw)

    def prepare_galaxy_input(self):
        '''
        writes the galaxy input file from a previously written template.
        Simply overwrites the filename.
        '''
        self.galaxy_inputs = []

        lines = open(self.galaxy_input).readlines()
        extra = ' '.join(lines[-3].split(' ')[1:])

        for i, sfr_file in enumerate(self.sfr_files):
            lines[-3] = ' '.join([self.sfr_files[i], extra])
            new_out = os.path.join(self.outfile_loc,
                                   os.path.split(self.galaxy_input)[1].replace('.dat', '%02i.dat' % i))
            with open(new_out, 'w') as f:
                f.write(''.join(lines))
            print 'wrote %s' % new_out
            self.galaxy_inputs.append(new_out)

    def write_results_table(self, trilgal_output):
        sgal = rsp.Galaxies.simgalaxy(trilegal_output, filter1='F606W', filter2='F814W')
        

    def vary_the_SFH(self, cmd_input_file, make_many_kw=None):
        make_many_kw = make_many_kw or {}
        self.prepare_trilegal_sfr(make_many_kw=make_many_kw)

        self.prepare_galaxy_input()

        agb_model = cmd_input_file.replace('cmd_input_', '').lower()
        target = os.path.split(self.galaxy_input)[1].replace('input_', '').replace('.dat', '').lower()

        trilegal_output = os.path.join(self.outfile_loc,
                                       'output_%s_%s' % (target, agb_model))
        for galaxy_input in self.galaxy_inputs:
            rsp.TrilegalUtils.run_trilegal(cmd_input_file, galaxy_input,
                                           trilegal_output, rmfiles=False,
                                           dry_run=True)
            self.write_results_table(trilegal_output)

def read_mettable():
    tab = os.path.join(table_src, 'IR_NAGBs.dat')
    dtype = [('fitstable', '|S46'),
             ('IR_TRGB', '<f8'),
             ('N_AGB', '<f8'),
             ('logOH', '<f8'),
             ('OHerr', '<f8'),
             ('Z', '<f8'),
             ('FeH', '<f8'),
             ('Ttype', '<f8')]
    table = np.genfromtxt(tab, dtype=dtype, delimiter=',')
    return table


def get_key_fromtable(ID, key):
    tab = read_mettable()
    names = list(tab['fitstable'])
    i, = [names.index(i) for i in names if ID in i]
    return tab[key][i]


def double_gaussian_optical_cmd_boxes(target=None):
    dg_box = {'DDO82': [0.6, 1.],
              'IC2574-SGS': [1., 1.8],
              'UGC-4305-1': [0.8, 1.6],
              'UGC-4305-2': [0.8, 1.6],
              'NGC4163': [0.5, 1.2],
              'UGC8508':  [0.5, 2.]}
    if target is None:
        return dg_box
    else:
        return dg_box[target]


def gaussian(x, p0, p1, p2):
    '''
    gaussian(arr,p): p[0] = norm, p[1] = mean, p[2]=sigma
    '''
    return p0 * np.exp( -1 * (x - p1) ** 2 / (2 * p2 ** 2))


def opt_cmd_contamination(target, band='opt'):
    import scipy
    from scipy import curve_fit
    gal = load_galaxy(target, band=band)
    points = np.column_stack((gal.Color, gal.Mag2))

    # stars in rgb box
    rgb_poly = get_opt_rgb_polygons(target)
    all_rgb_stars, = np.nonzero(nxutils.points_inside_poly(points, rgb_poly))
    
    # stars in the rheb and rgb region
    dg_box = double_gaussian_optical_cmd_boxes(target=target)
    col_bins = np.arange(dg_box[0], dg_box[1], 0.05)
    mag_bins = np.arange(gal.Trgb, gal.Trgb + 1.5, 0.15)
    
    # Define rheb+rgb boxes
    rheb_rgb_starss = []
    rgb_starss = []
    fig, ax = plt.subplots()
    for u, l in zip(mag_bins, np.roll(mag_bins,-1))[:-1]:
        verts = np.array([[dg_box[0], l],
                          [dg_box[0], u],
                          [dg_box[1], u],
                          [dg_box[1], l],
                          [dg_box[0], l]])
        rheb_rgb_stars, = np.nonzero(nxutils.points_inside_poly(points, verts))
        rheb_rgb_starss.append(rheb_rgb_stars)
        rgb_stars = list(set(all_rgb_stars) & set(rheb_rgb_stars))
        rgb_starss.append(rgb_stars)
    
    for i in range(len(rheb_rgb_starss)):
        # stars in rheb+rgb and not in rgb box
        rgb_stars = rgb_starss[i]
        rheb_stars = list(set(rheb_rgb_starss[i]) - set(rgb_stars))
        ax.plot(gal.Color[rheb_stars], gal.Mag2[rheb_stars], '.')
        ax.plot(gal.Color[rgb_stars], gal.Mag2[rgb_stars], '.')

        mag_hist = np.histogram(gal.Color[rheb_stars], bins=col_bins)[0]
        # fit gaussian to rheb_stars
        try:
            parameters, covariance = curve_fit(gaussian, col_bins[1:],
                                               mag_hist, p0=[np.max(mag_hist), 1, 0.5])
        except RuntimeError:
            print 'can not fit %i' % i
        big_color_array = np.arange(col_bins[0], col_bins[-1], 0.001)
        fit_gauss = gaussian(big_color_array, *parameters)
        fig, ax = plt.subplots()
        ax.plot(big_color_array, fit_gauss)
        ax.plot(col_bins[1:], mag_hist, ls='steps--')
        rgb_mag_hist = np.histogram(gal.Color[rgb_stars], bins=col_bins)[0]
        ax.plot(col_bins[1:], rgb_mag_hist, ls='steps--')
        ax.set_title(mag_bins[i])
        
        # see how many stars from rheb are in rgb box
        rgb_blue_edge = col_bins[np.max(np.nonzero(mag_hist > 1)) + 1]
        rheb_in_rgb, err = scipy.integrate.quad(gaussian, rgb_blue_edge, np.inf,
                                                args=(parameters[0], parameters[1],
                                                      parameters[2]))
        
        
        
        ax.set_title('%.2f, %.2f' % (mag_bins[i], rheb_in_rgb))
        
        
        # now do the same with rgb_stars and see how many are redder than the rgb
        # blue edge
    #6. take percentage of 5, adjust the edge of 1 closer if necessary.
    
    
    
    '''
    for u, l in zip(mag_bins, np.roll(mag_bins,-1))[:-1]:
        print u
        verts = np.array([[dg_box[0], l],
                          [dg_box[0], u],
                          [dg_box[1], u],
                          [dg_box[1], l],
                          [dg_box[0], l]])
        ax.plot(verts[:, 0], verts[:, 1])
        rstars_in_agb, agb_in_rstars, poisson_noise, nstars, color_sep = \
        gal.double_gaussian_contamination(verts, diag_plot=True,   
                                                   dcol=0.05, Color=gal.Color,
                                                   Mag2=gal.Mag2)
    '''
    '''
    targets = load_targets('paper1')
    gals = rsp.Galaxies.galaxies([load_galaxy(t, band=band) for t in targets])
    g555 = rsp.Galaxies.galaxies(gals.select_on_key('filter1', 'F555W'))
    g606 = rsp.Galaxies.galaxies(gals.select_on_key('filter1', 'F606W'))
    g475 = rsp.Galaxies.galaxies(gals.select_on_key('filter1', 'F475W'))
    
    for g in [g555, g606, g475]:
        g.squish('Color', 'Mag2', 'Trgb')
        offset = gals.Trgbs - np.mean(gals.Trgbs)
        Mag2 = np.concatenate([gg.Mag2 - offset[i] for i, gg in enumerate(g.galaxies)])
        Color = g.Colors
        
    
    '''

    
def cmd_contamination(targets, band='ir', dcol=0.02, dmag=0.2, thresh=5):
    '''
    , verts_file_fmt=None,
    diag_plot=True, isoch_base=None,
    isoc_search_fmt=None, out_file_fmt=None,
    out_fig_fmt=None, M11=False, galaxy_table=None,
    rdcolor=None, rdmag=None
    '''
    targets = load_targets(targets)
    # load galaxies class
    gals = rsp.Galaxies.galaxies([load_galaxy(t, band=band) for t in targets])
    # combine all galaxies
    gals.squish('Color', 'Mag2', 'Trgb')    
    mean_trgb = np.mean(gals.Trgbs)
    offset = gals.Trgbs - np.mean(gals.Trgbs)
    Mag2 = np.concatenate([g.Mag2 - offset[i] for i, g in enumerate(gals.galaxies)])
    Color = gals.Colors
    gals.Mag2o = Mag2
    verts = agb_rheb_separation()
    fig, ax = gals.galaxies[0].plot_cmd(gals.Colors, gals.Mag2o,
                                        threshold=75, levels=5,
                                        hist_bin_res=.1)
    # rough outline of the RGB/AGB area. Exclude MS (~0.0)
    ms_cut = 0.6
    # minimum size of dmag bin spacing

    #rstars_in_agb, agb_in_rstars, poisson_noise, nstars, color_sep = \
    #    gals.galaxies[0].double_gaussian_contamination(verts, diag_plot=True,   
    #                                                   dcol=0.02, Color=Color,
    #                                                   Mag2=Mag2)
    points = np.column_stack((Color, Mag2))
    agb_stars, = np.nonzero(nxutils.points_inside_poly(points, verts))
    rgb_rheb_stars = list(set(np.nonzero((Color > ms_cut) &
                                          (Mag2 > np.min(verts[:,1])) &
                                          (Mag2 < np.max(verts[:,1])))[0]) - \
                          set(agb_stars))
    out = open('tpagb_contamination.dat', 'w')
    mag_bins = np.arange(-10, mean_trgb + dmag,  dmag)
    for u, l in zip(mag_bins, np.roll(mag_bins,-1))[:-1]:
        cbinsl = np.arange(ms_cut, np.max(Color[rgb_rheb_stars]), dcol)
        mag_left, = np.nonzero((Mag2[rgb_rheb_stars] > u) &
                               (Mag2[rgb_rheb_stars] < l))

        mag_right, = np.nonzero((Mag2[agb_stars] > u) &
                                (Mag2[agb_stars] < l))
        cbinsr = np.arange(np.min(Color[agb_stars][mag_right]),
                           np.max(Color[agb_stars][mag_right]), dcol)
        
        ax.plot(Color[rgb_rheb_stars][mag_left], Mag2[rgb_rheb_stars][mag_left],'o')
        ax.plot(Color[agb_stars][mag_right], Mag2[agb_stars][mag_right], 'o')
        mag_hist_left, _ = np.histogram(Color[rgb_rheb_stars][mag_left],
                                        bins=cbinsl)
        mag_hist_right, _ = np.histogram(Color[agb_stars][mag_right],
                                         bins=cbinsr)
        linr, rinl = HeB_contamination(cbinsl, mag_hist_left, cbinsr, mag_hist_right,
                          diag_plot=True, ax=ax, mag=np.mean([u, l]))
        out.write('%.3f %.3f %.3f \n' % (u, linr/float(np.sum(mag_hist_right)), rinl/float(np.sum(mag_hist_left))))
    out.close()


def HeB_contamination(color_left, mag_hist_left, color_right, mag_hist_right,
                      diag_plot=False, ax=None, mag=0):
    '''
    This function fits a double gaussian to a color histogram of stars
    within the <maglimits> and <colorlimits> (tuples).

    It then finds the intersection of the two gaussians, and the fraction
    of each integrated gaussian that crosses over the intersection color
    line.
    '''
    try:
        mpfit
    except NameError:
        from mpfit import mpfit
        
    try:
        integrate
    except NameError:
        from scipy import integrate
    # the indices of the stars within the MS/BHeB regions
    import functions
    # uniform errors
    erra = np.zeros(len(color_left[:1])) + 1.
    errb = np.zeros(len(color_right[:1])) + 1.

    # set up inputs
    hist_ina = {'x': color_left[1:], 'y': mag_hist_left, 'err': erra}
    hist_inb = {'x': color_right[1:], 'y': mag_hist_right, 'err': errb}
    
    # set up initial parameters:
    # norm = max(hist),
    # mean set to be half mean, and 3/2 mean,
    # sigma set to be same as dcol spacing...
    p0a = [np.nanmax(mag_hist_left), np.mean(color_left) / 2,
           np.diff(color_left)[0]]
    p0b = [np.nanmax(mag_hist_right), np.mean(color_right) / 2,
           np.diff(color_right)[0]]

    #mp_dg = mpfit(functions.mp_double_gauss, p0, functkw=hist_in, quiet=True)
    mp_dga = mpfit(functions.mp_gaussian, p0a, functkw=hist_ina, quiet=True)
    mp_dgb = mpfit(functions.mp_gaussian, p0b, functkw=hist_inb, quiet=True)
    #
    #if mp_dg.covar is None:
    #    return 0., 0., 0., 0.
    #perc_err = (np.array(mp_dg.perror) - np.array(mp_dg.params)) / \
    #            np.array(mp_dg.params)
    #if np.sum([p ** 2 for p in perc_err]) > 10.:
    #    return 0., 0., 0., 0.
    # take fit params and apply to guassians on an arb color scale
    color_array = np.linspace(-10, 10, 1000)
    #g_p1 = mp_dg.params[0: 3]
    #g_p2 = mp_dg.params[3:]
    #gauss1 = functions.gaussian(color_array, g_p1)
    #gauss2 = functions.gaussian(color_array, g_p2)
    gauss1 = functions.gaussian(color_array, mp_dga.params)
    gauss2 = functions.gaussian(color_array, mp_dgb.params)
    # color separatrion is the intersection of the two gaussians..
    # min_locs = rsp.math_utils.find_peaks(gauss1 + gauss2)['minima_locations']
    auto_color_sep = [] #color_array[min_locs]
    if len(auto_color_sep) == 0:
        auto_color_sep = np.mean([np.max(color_left), np.min(color_right)])
    print np.max(color_left), auto_color_sep, np.min(color_right)
    # find contamination past the color sep...
    #g12_Integral = integrate.quad(functions.double_gaussian, -np.inf, np.inf, mp_dg.params)
    #try:
    #    norm =  float(len(all_inds)) / g12_Integral[0] 
    #except ZeroDivisionError:
    #    norm = 0.
    g1_Integral = integrate.quad(functions.gaussian, -np.inf, np.inf,
                                 mp_dga.params)
    g2_Integral = integrate.quad(functions.gaussian, -np.inf, np.inf,
                                 mp_dgb.params)
    norm = g1_Integral[0] + g2_Integral[0]
    if norm == 0:
        norm = 1.
    g1_Int_colsep = integrate.quad(functions.gaussian, -np.inf, auto_color_sep,
                                   mp_dga.params)
    g2_Int_colsep = integrate.quad(functions.gaussian, auto_color_sep, np.inf,
                                   mp_dgb.params)

    left_in_right = (g1_Integral[0] - g1_Int_colsep[0])
    right_in_left = (g2_Integral[0] - g2_Int_colsep[0])
    '''
    
    try:
        left_in_right = g1_Int_colsep[0] / g1_Integral[0]
    
        left_in_right = 0.

    try:
        right_in_left = g2_Int_colsep[0] / g2_Integral[0]
    except ZeroDivisionError:
        right_in_left = 0.
    '''
    # diagnostic
    #print color_sep
    if diag_plot is True:
        if ax is None:
            fig, ax = plt.subplots()
        ax.plot(color_left[1:], mag_hist_left/np.sum(mag_hist_left) + mag, ls='steps', lw=2, color='navy')
        ax.plot(color_right[1:], mag_hist_right/np.sum(mag_hist_left) + mag, ls='steps', lw=2, color='darkred')
        #ax1.plot(color_array,
        #        functions.double_gaussian(color_array, mp_dg.params))
        ax.plot(color_array, gauss1//np.sum(gauss1)+gauss2/np.sum(gauss2) + mag,
                color='purple', lw=2)
        ax.plot(color_array, gauss1/np.sum(gauss1) + mag, color='blue', lw=2)
        ax.plot(color_array, gauss2/np.sum(gauss2) + mag, color='red', lw=2)
        #ax1.set_ylim((0, 100))
        #ax.set_xlim(-1, 2.5)
        #ax1.set_xlabel('$%s-%s$' % (gal.filter1, gal.filter2), fontsize=20)
        #ax.set_ylabel('$\#$', fontsize=20)
        #ax1.set_title('$%s$' % gal.target)
        ax.vlines(auto_color_sep, *ax.get_ylim())
        ax.text(-0.5, mag, 'left in right: %i' % left_in_right)
        ax.text(2, mag, 'right in left: %i' % right_in_left, ha='right')
        #fig1.savefig('heb_contamination_%s_%s_%s_mag2_%.2f.png' % (gal.filter1, gal.filter2, gal.target,np.mean(np.array(all_verts)[:, 1])))
        #plt.close()
    return left_in_right, right_in_left


def get_imf(target):
    filename = os.path.join(table_src, 'best_fits_from_match_runs.dat')
    dtype = [('ID', '|S5'), ('Galaxy', '|S14'), ('filter1', '|S5'),
             ('filter2', '|S5'), ('Av', '<f8'), ('IMF', '<f8'),
             ('dmod', '<f8'), ('dlogZ', '<f8')]
    data = np.genfromtxt(filename, dtype=dtype)
    imf, = data['IMF'][np.nonzero(data['Galaxy']==target)]
    if imf == 1.30:
        tab_imf = 'tab_imf/imf_salpeter_match.dat'        
    elif imf == 1.35:
        tab_imf = 'tab_imf/imf_salpeter.dat'
    else:
        print 'imf problem.'
    return tab_imf


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


def plot_opt_hess():
    targets = targets_paper1()
    gals = rsp.Galaxies.galaxies([load_galaxy(t, band='opt') for t in targets])

    g606 = gals.select_on_key('filter1','F606W')
    '''
    poly_dict = {'DDO82': np.array([[ 0.6232878 , -2.96601422],
                                 [ 0.74828178, -3.41942444],
                                 [ 0.82520115, -3.82245574],
                                 [ 0.90212052, -4.12472922],
                                 [ 1.62323962, -4.05755733],
                                 [ 1.52709041, -3.63773306],
                                 [ 1.45978596, -3.38583849],
                                 [ 1.38286659, -3.11715096],
                                 [ 1.30594722, -2.86525639],
                                 [ 1.23864277, -2.59656886],
                                 [ 0.50790874, -2.56298292]]),
                 'NGC4163': array([[ 0.59977447, -2.57185878],
                                   [ 0.746842  , -2.99357265],
                                   [ 0.89390953, -3.52372723],
                                   [ 0.94035191, -4.0177349 ],
                                   [ 1.28866974, -4.04183284],
                                   [ 1.18030419, -3.60807   ],
                                   [ 1.14160221, -3.21045407],
                                   [ 1.09515983, -2.62005465]])}
    '''

    g606s = rsp.Galaxies.galaxies(g606)
    g606s.squish('Color', 'Mag2', 'Trgb')

    mean_trgb = np.mean(g606s.Trgbs)
    offset = g606s.Trgbs - np.mean(g606s.Trgbs)
    g606s.Mag2o = np.concatenate([g.Mag2 - offset[i] for i, g in enumerate(g606s.galaxies)])

    g606hess =  rsp.astronomy_utils.hess(g606s.Colors, g606s.Mag2o, 0.1, cbinsize=0.05)


    extent = [g606hess[0][0], g606hess[0][-1], g606hess[1][-1], g606hess[1][0]]
    imshow_kw={'norm': LogNorm(vmin=None, vmax=g606hess[2].max()),
               'cmap': plt.cm.gray_r, 'interpolation': 'nearest',
               'aspect': 'equal', 'extent': extent}

    fig, ax = plt.subplots()
    ax.autoscale(False)
    ax.set_xlim(-1, 3)
    ax.set_ylim(0.5, -7)
    ax.imshow(g606hess[2], **imshow_kw)
    ax.set_aspect(1./2.)
    ax.set_xlabel('$F606W-F814W$', fontsize=20)
    ax.set_ylabel('$F814W$', fontsize=20)
    ax.hlines(mean_trgb, *ax.get_xlim(), lw=2, color='red', zorder=100)
    #ax.hlines(mean_trgb+1.5, *ax.get_xlim(), lw=2, color='red', zorder=100)
    #for k, v in poly_dict.items():
    #    ax.plot(v[:,0], v[:,1])
    #plt.colorbar(cs)
    plt.savefig('opt_cmd_f606w.png', dpi=300)


def multi_galaxy_hess(targets=None, split_by_color=False, ax=None, imshow_kw={},
                      make_hess=False, band='ir', fits_src='default'):
    if ax is None:
        fig, ax = plt.subplots()

    targets = load_targets(targets)
    # load galaxies class
    gals = rsp.Galaxies.galaxies([load_galaxy(t, band=band, fits_src=fits_src) for t in targets])
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
        cmin, cmax, cbinsize = -0.25, 1.25, 0.01
        if make_hess is True:
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
                                                threshold=75, levels=5,
                                                hist_bin_res=.1)
            #fig.colorbar(gals.galaxies[0].cs)
        verts = agb_rheb_separation()
        ax.plot(verts[:, 0], verts[:, 1], lw=2, color='navy')
        ax.hlines(mean_trgb, *ax.get_xlim(), lw=2, color='red', zorder=100)
        ax.set_xlim(cmin, cmax)

    ax.set_xlabel('$%s-%s$' % (gals.filter1s[0], gals.filter2s[0]), fontsize=20)
    ax.set_ylabel('$%s$' % gals.filter2s[0], fontsize=20)
    ax.tick_params(labelsize=16)
    return fig, ax, gals


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
    vert_file = research_path + 'code/TPAGB-calib/agb_rheb_sep.dat'
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


def load_galaxy(ID, band='ir', fits_src='default'):
    '''
    '''
    if fits_src == 'default':
        fits_src = os.path.join(snap_src, 'data', 'galaxies')

    if band is None:
        fits_src = os.path.join(snap_src, 'data', 'opt_ir_matched_v2')

    fitsname = rsp.fileIO.get_files(fits_src, '*%s*fits' % ID)
    hla = False
    if len(fitsname) == 0:
        fitsname = rsp.fileIO.get_files(fits_src, '*%s*fits' % ID.lower())
        hla = True

    if band == 'opt':
        filetype = 'fitstable'
    else:
        filetype = 'agbsnap'

    fitsname = ir_or_opt_file(fitsname, band=band)

    gal_kw = {'hla': hla, 'photsys': 'wfc3snap', 'angst': True,
              'filetype': filetype, 'band': band}

    gal_kw['z'] = get_key_fromtable(ID, 'Z')

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


def get_sfr_file(ID, sfr_dir='default'):
    if sfr_dir == 'default':
        sfr_dir = os.path.join(snap_src, 'data', 'sfh')
    return rsp.fileIO.get_files(sfr_dir, '*%s*' % ID)[0]


def compare_metallicities():
    IDs = all_targets()
    zs = match_metallicities(IDs)
    for ID in IDs:
        zs[ID]['zmeas'] = get_key_fromtable(ID, 'Z')

    for id, zdict in sorted(zs.items(), key=lambda(k, v): (v['avez'], k)):
        print '%s %.4f %.4f' % (id, zdict['avez'], zdict['zmeas'])


def setup_trilegal(gal, model, object_mass=5e9, sfr_dir='default',
                   sfr_file=None, outfile_loc='default'):
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
    if outfile_loc == 'default':
        outfile_loc = snap_src
    agb_model = model.replace('cmd_input_', '').replace('.dat', '').lower()

    object_dist = 10 ** ((5 + gal.dmod) / 5)
    if sfr_file is None:
        sfr_file = get_sfr_file(gal.target, sfr_dir=sfr_dir)

    gal_dict_inp = {'photsys': gal.photsys,
                    'filter1': gal.filter2,
                    'object_mass': object_mass,
                    'object_sfr_file': sfr_file}

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

    galaxy_input = os.path.join(outfile_loc, 'input', 'input_%s.dat' % gal.target)
    gal_inp.write_params(galaxy_input, rsp.TrilegalUtils.galaxy_input_fmt())
    trilegal_output = os.path.join(outfile_loc, 'output', 'output_%s_%s.dat' %
                                   (gal.target, agb_model))

    return galaxy_input, trilegal_output, galinp_kw, gal_inp


def setup_data_normalization(gal, filt1, filt2, band=None, leo_method=False,
                             use_opt_rgb=False):
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
        # the 2 std around color mean
        cmean = np.mean(color[rgb_norm])
        cstd = np.std(color[rgb_norm]) * 2.
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
    if use_opt_rgb is True:
        ndata_stars = get_opt_nrgb(gal.target)
        gal.rgb_norm_inds = ndata_stars
    else:
        gal.rgb_norm_inds = list(set(rgb_norm) & set(inds))
        ndata_stars = len(gal.rgb_norm_inds)

    gal.norm_verts = verts
    logger.debug('norm verts:')
    logger.debug(verts)
    logger.info('data stars in normalization range: %i' % ndata_stars)
    return verts, ndata_stars


def run_make_normalized_simulation(targets, models, band='ir', sfr_dir='default'):
    targets = load_targets(targets)

    if type(models) == str:
        models = [models]

    for target in targets:
        gal = load_galaxy(target, band)
        print target
        for model in models:
            print model
            make_normalized_simulation(gal, model, gal.filter1, gal.filter2,
                                       object_mass=load_sim_masses(gal.target),
                                       band=band, sfr_dir=sfr_dir)


def make_normalized_simulation(gal, model, filt1, filt2, photsys='wfc3snap',
                               over_write=False, object_mass=5e6,
                               maglims='trgb', band='ir', offsets=(1.5, 0.5),
                               count_offset=0., mass_inc_fact=5,
                               run_trilegal=True, norm_threshold=0.75,
                               leo_method=False, spread_outfile=None,
                               leo_norm=False, leo_ast=False,
                               spread_outfile2=None, trilegal_output=None,
                               norm_fname=None, use_opt_rgb=False,
                               sfr_dir='default'):
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
        galaxy_input, trilegal_output, galinp_kw, gal_inp = \
            setup_trilegal(gal, model, object_mass=object_mass, sfr_dir=sfr_dir)
        if band == 'opt':
            trilegal_output = trilegal_output.replace('.dat', '_opt.dat')

        cmd_input = os.path.join(research_path + 'padova_apps/cmd_inputfiles/',
                                 model)

    # might not be needed, but better to load them outside the while loop.
    fake_files = get_fake_files(gal.target, band=band)

    # set maglim if not already set.
    if maglims == 'trgb':
        gal.maglims = (trgb + offsets[0], trgb + offsets[1])
        logger.info('normalization mag: %f below the trgb to %f above' % (offsets[0], offsets[1]))
    else:
        logger.info('normalization mag: %f to %f' % (offsets[0], offsets[1]))

    verts, ndata_stars = setup_data_normalization(gal, filt1, filt2,
                                                  leo_method=leo_norm,
                                                  use_opt_rgb=use_opt_rgb)

    # initializations
    go = 0
    normalization = 1e9

    while normalization > norm_threshold:
        if go > 0:
            # if we've been through this already, increase mass and try again.
            object_mass *= mass_inc_fact
            galinp_kw['object_mass'] = object_mass
            gal_inp.add_params(galinp_kw)
            gal_inp.write_params(galaxy_input, rsp.TrilegalUtils.galaxy_input_fmt())
        go += 1

        if run_trilegal is True:
            # run trilegal
            logger.debug('Trying %s %s, Mass=%g, Attempt %i' %
                         (gal.target, model, object_mass, go))

            rsp.TrilegalUtils.run_trilegal(cmd_input, galaxy_input,
                                           trilegal_output, rmfiles=False)

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
                               maglims=None, offsets=None, leo_norm=False,
                               sfr_dir='default'):

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
                                      trilegal_output=sim_file, sfr_dir=sfr_dir)

    return gal, sgal


def get_opt_nrgb(target):
    offset_dict = {'SCL-DE1': 1554,
                   #'NGC2403-HALO-6': -0.155,
                   #'NGC7793-HALO-6': -0.307,
                   #'UGC-04459': 0.063,
                   'UGC-4305-1': 6242,
                   'UGC-4305-2': 7814,
                   #'UGC-5139': 0.105,
                   'DDO78': 6542,
                   'DDO82': 23796,
                   #'KDG73': 0.157,
                   #'KKH37': 0.085,
                   #'NGC0300-WIDE1': -0.212,
                   #'NGC2403-DEEP': -0.089,
                   #'NGC2976-DEEP': -0.220,
                   #'NGC3741': 0.143,
                   'NGC4163': 10747,
                   'UGC8508': 3047,
                   #'UGCA292': 0.394,
                   #'NGC3077-PHOENIX': -0.297,
                   'IC2574-SGS': 15823,
                   #'HS117': 0.092,
                   'DDO71': 3665}
    return float(offset_dict[target])


def get_opt_rgb_polygons(target):
    '''
    I ran get_nrgb_from_optical with default args and pasted the ipython output here
    # target filter1 filter2 NRGB_box NRGB_by_eye polygon
    DDO78 F475W F814W 3842 3842 
    DDO71 F606W F814W 2248 2248 
    SCL-DE1 F606W F814W 944 944 
    DDO82 F606W F814W 14314 14331     
    IC2574-SGS F555W F814W 10485 8903 
    UGC-4305-1 F555W F814W 5337 3654  
    UGC-4305-2 F555W F814W 6555 4402  
    NGC4163 F606W F814W 7010 6136     
    UGC8508 F475W F814W 2411 2078     
    '''
    poly_dict = {'DDO78': 3842,
                 'DDO71': 2248,
                 'SCL-DE1': 944,
                 'DDO82': np.array([(0.62328780081964297, -2.9660142201464623), (0.74828177882779023, -3.4194244355894625), (0.82520114990972671, -3.8224557382054627), (0.90212052099166362, -4.1247292151674628), (1.6232396248848207, -4.0575573313981295), (1.5270904110323995, -3.6377330578397959), (1.4597859613357049, -3.3858384937047958), (1.3828665902537685, -3.1171509586274624), (1.3059472191718315, -2.8652563944924623), (1.238642769475137, -2.5965688594151288), (0.50790874419673759, -2.5629829175304621)]),
                 'IC2574-SGS': np.array([(0.88085277986008848, -2.594209765189702), (0.97055019845183477, -2.8857664759318862), (1.1179102432811323, -3.2502123643596157), (1.2012007034020398, -3.5539172713827236), (1.2844911635229472, -3.8333257858439822), (1.354967706702177, -4.0034005337769232), (1.3934094575272109, -4.0641415151815448), (1.9700357199027243, -4.0398451226196963), (1.8290826335442656, -3.6025100565064205), (1.6625017133024507, -3.2137677755168426), (1.6368738794190945, -2.9586556536174315), (1.5728042947107044, -2.6185061577515505), (1.5599903777690258, -2.5820615689087782)]),
                 'UGC-4305-1': np.array([(0.96940535804821337, -2.6403426463456499), (1.0757322587834022, -3.0381712503347198), (1.182059159518591, -3.4239444420816976), (1.2508589188178307, -3.6409418624393717), (1.3321677252623871, -3.8217730460707675), (1.4134765317069435, -4.0628812909126282), (1.5385670031601064, -4.0990475276389073), (1.7824934224937756, -4.0990475276389073), (1.6636574746132702, -3.7735513971023953), (1.6073667624593462, -3.4721660910500698), (1.5510760503054231, -3.1105037237872781), (1.476021767433525, -2.7729521810086735), (1.4572581967155509, -2.5559547606509985)]),
                 'UGC-4305-2': np.array([(0.9928229629548313, -2.6090931485301443), (1.0590397476993938, -2.8305905745508753), (1.1071974093318029, -3.1013096507984343), (1.1613747786682622, -3.2858908391490429), (1.2336112711168759, -3.6796640409636749), (1.2817689327492849, -3.8888560544276984), (1.3359463020857452, -4.0611318302215995), (1.4262419176465113, -4.098048067891721), (1.6008134410639934, -4.122658893005136), (1.7091681797369138, -4.1103534804484285), (1.7814046721855274, -4.073437242778307), (1.7212075951450161, -3.6796640409636749), (1.6369316872883002, -3.1751421261386783), (1.5707149025437386, -2.8182851619941678), (1.5285769486153806, -2.6213985610868518)]),
                 'NGC4163': np.array([(0.59977446601686069, -2.5718587825132402), (0.7468419971858129, -2.9935726500356323), (0.89390952835476511, -3.523727226349497), (0.94035190661864476, -4.0177348997328703), (1.2886697435977403, -4.0418328350198642), (1.1803041943153563, -3.6080699998539747), (1.1416022124287899, -3.2104540676185764), (1.0951598341649103, -2.620054653087228)]),
                 'UGC8508': np.array([(1.5332306001668927, -2.5908838137984285), (1.9516553658152, -3.9626182764768609), (2.0258920177850603, -4.0681363120675096), (2.4510655699760822, -4.0564120858907708), (2.5590461546595167, -3.9977909550070772), (2.2958434794936462, -3.4936492294073114), (2.093379883212207, -2.9777832776308069), (1.9246602196443412, -2.5557111352682123)])
                 }
    return poly_dict[target]


def get_ir_rgb_polygons(gal, leo_method=False):
    if leo_method is False:
        # find normalization region
        # (when hand picking cmd_regions, RGB includes AGB!)
        irgb, = np.nonzero(gal.color > 0.4)
        logger.info('normalization color: using the stdev around the mean of stars in mag range using color cut of 0.4')
        # the stars in the data marked as rgb between mag lims.
        rgb_norm = rsp.math_utils.between(gal.mag2, gal.maglims[0], gal.maglims[1],
                                          inds=irgb)
        # the 2 std around color mean
        cmean = np.mean(gal.color[rgb_norm])
        cstd = np.std(gal.color[rgb_norm]) * 2.
        col_min = cmean - cstd
        col_max = cmean + cstd
    else:
        # Leo's method is to use all the stars between the mag lims.
        logger.info('normalization color: using all stars in mag range')
        rgb_norm = rsp.math_utils.between(gal.mag2, gal.maglims[0], gal.maglims[1])
        col_min = np.min(gal.color)
        col_max = np.max(gal.color)

    verts = np.array([[col_min, gal.maglims[0]],
                      [col_min, gal.maglims[1]],
                      [col_max, gal.maglims[1]],
                      [col_max, gal.maglims[0]],
                      [col_min, gal.maglims[0]]])

    return verts


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
            'NGC300-WIDE1',
            'NGC3077-PHOENIX',
            'NGC2403-DEEP',
            'NGC2403-HALO-6',
            'UGC-5139',
            'DDO82',
            'IC2574-SGS',
            'UGC4305-1',
            'UGC4305-2',
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
    #return ['UGC-4305-1', 'UGC-4305-2', 'NGC4163', 'UGC8508']
    return np.concatenate([gi10_overlap(), targets_z002()])


def targets_leo_norm_ok():
    '''
    targets that don't have MS to cut out, might as well use the whole
    RGB
    '''
    return np.concatenate([gi10_overlap(), ['DDO82', 'NGC4163', 'UGC8508']])


def add_opt_asts_to_ir_asts(targets, models, ast_file=None,
                            do_ast_ir=False):
    targets = load_targets(targets)
    for target in targets:
        opt_fake = get_fake_files(target, band='opt')
        print opt_fake
        opt_gal = load_galaxy(target, band='opt')
        print opt_gal
        ir_gal = load_galaxy(target)
        for model in models:
            if ast_file is not None:
                ir_ast_file = ast_file
            else:
                ir_ast_file = load_ast_file(ir_gal, model)
            out_file = '_'.join(np.concatenate([ir_ast_file.split('_')[:3],
                                opt_gal.filters, ir_ast_file.split('_')[3:]]))

            sgal = rsp.Galaxies.simgalaxy(ir_ast_file, filter1=opt_gal.filter1,
                                          filter2=opt_gal.filter2)
            print 'four filter ast: %s' % out_file
            rsp.Galaxies.ast_correct_trilegal_sim(sgal, fake_file=opt_fake,
                                                  outfile=out_file)

    
def load_ast_file(gal, model, extra=''):
    if not 'F1' in gal.filter2:
        extra = '_opt'
    model_short = model.replace('cmd_input_','').lower().replace('.dat','%s.dat' % extra)
    search_str = '*'.join(('ast_%s' % gal.filter1, gal.filter2 + '_output',
                           gal.target, model_short))
    try:
        sim_file, = rsp.fileIO.get_files(os.path.join(snap_src, 'output'),
                                         search_str)
    except ValueError:
        try:
            search_str = '*'.join(('ast', gal.filter1, gal.filter2 + '_output',
                               gal.target, model_short))
            sim_file, = rsp.fileIO.get_files(os.path.join(snap_src, 'output'),
                                             search_str)
            print 'using four filter.'
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


def model_ratio_once(target, model, ast_file_loc, rgb_verts, ndata_rgb, agb_verts,
                     opt_agb_verts, nar_ratio_data, nar_ratio_data_opt, ir_gal,
                     opt_gal, find_file, figname, ndata_agb, plot_LF_kw,
                     ir_data_pts, gal_ir_hist, gal_ir_bins, gal_opt_hist,
                     gal_opt_bins, data_pts, ndata_agb_opt):
    logger.info('%s: %s' % (target, model))
    model_short = model.replace('cmd_input_', '').lower()
    # load four filter ast file
    find_file['model'] = model_short 
    search_str = 'ast_%(filter1)s_%(filter2)s_%(filter3)s_%(filter4)s_output_%(target)s_%(model)s' % find_file
    sim_file, = rsp.fileIO.get_files(ast_file_loc, search_str)

    sgal = rsp.Galaxies.simgalaxy(sim_file, filter1='F110W', filter2='F160W')
    sgal.load_ast_corrections()
    logger.info('%s: %s loaded asts.' % (target, model))
    # so i dunno why, but filter3 is 814, filter4 is the other guy.
    opt_sim_pts = np.column_stack((sgal.ast_mag4 - sgal.ast_mag3, sgal.ast_mag3))
    ir_sim_pts = np.column_stack((sgal.ast_mag1 - sgal.ast_mag2, sgal.ast_mag2))
    
    # make normalization from the full simulation
    rgb_in_sim, = np.nonzero(nxutils.points_inside_poly(opt_sim_pts,
                                                        rgb_verts))
    nsim_rgb = float(len(rgb_in_sim))
    normalization = ndata_rgb / nsim_rgb

    nmodel = np.array([])
    nmodel_opt = np.array([])
    # do lots of random sampling.
    logger.info('%s: %s going into the mc loop.' % (target, model))
    for i in range(10001):
        # random sample the data distribution
        rands = np.random.random(len(sgal.ast_mag4))
        ind, = np.nonzero(rands < normalization)
        if i % 5000 == 0: 
            logger.info('%s: %s mc loop... %i' % (target, model, i))
        norm_sim_pts_opt = opt_sim_pts[ind]
        norm_sim_pts_ir = ir_sim_pts[ind]

        # number of agb and rgb stars in the scaled simulation
        sim_agb_norm, = np.nonzero(nxutils.points_inside_poly(norm_sim_pts_ir,
                                                              agb_verts))
        sim_rgb_norm, = np.nonzero(nxutils.points_inside_poly(norm_sim_pts_opt,
                                                              rgb_verts))
        sim_agb_norm_opt, = np.nonzero(nxutils.points_inside_poly(norm_sim_pts_opt,
                                                                  opt_agb_verts))

        # simulation's nagb/nrgb ratio
        nsim_agb_opt = float(len(sim_agb_norm_opt))
        nsim_agb = float(len(sim_agb_norm))
        nsim_rgb = float(len(sim_rgb_norm))
        nar_ratio_sim = nsim_agb / nsim_rgb

        nar_ratio_sim_opt = nsim_agb_opt / nsim_rgb
        nmodel = np.append(nmodel, nar_ratio_sim)
        nmodel_opt = np.append(nmodel_opt, nar_ratio_sim_opt)
    logger.info('%s: %s finished mc loop.' % (target, model))
    mc_pct_diff = (np.mean(nmodel) - nar_ratio_data) / nar_ratio_data
    mc_pct_diff_opt = (np.mean(nmodel_opt) - nar_ratio_data_opt) / nar_ratio_data_opt

    sgal_ir_hist, sbins_ir = hist_it_up(norm_sim_pts_ir[:, 1])
    sgal_opt_hist, sbins_opt = hist_it_up(norm_sim_pts_opt[:, 1])

    # use mean nmodels for plots with lines.

    # plot the IR LF
    model_title = translate_model_name(model)
    plot_LF_kw['model_title'] = model_title
    logger.info('%s: %s plotting IR' % (target, model))
    fig, axs, _ = ir_gal.plot_LF(ir_data_pts[:,0], ir_data_pts[:,1],
                 norm_sim_pts_ir[:, 0], norm_sim_pts_ir[:, 1],
                 ir_gal.filter1, ir_gal.filter2,
                 gal_hist=gal_ir_hist, bins=gal_ir_bins,
                 sgal_hist=sgal_ir_hist, sbins=sbins_ir,
                 **plot_LF_kw)
    logger.info('%s: %s IR plotted' % (target, model))
    add_lines_LF(fig, axs, ir_gal, [ir_gal.trgb, ir_gal.trgb+2],
                 [ir_gal.model_plt_color, ir_gal.data_plt_color],
                 [ndata_agb, nsim_agb, ndata_rgb, nsim_rgb])  
    
    fig_name_kw = {'target': opt_gal.target,
                   'model_name': model.replace('cmd_input_', '').replace('.dat', ''),
                   'filter1': ir_gal.filter1,
                   'filter2': ir_gal.filter2}
    fig.savefig(figname % fig_name_kw, dpi=300,
                bbox_to_inches='tight')
    fig.close()
    
    if target in gi10_overlap():
        # if it's part of G10 plot the OPT LF
        logger.info('%s: %s plotting OPT' % (target, model))
        fig, axs, _ = opt_gal.plot_LF(data_pts[:,0], data_pts[:,1],
                       norm_sim_pts_opt[:, 0], norm_sim_pts_opt[:, 1],
                       opt_gal.filter1, opt_gal.filter2,
                       gal_hist=gal_opt_hist, bins=gal_opt_bins,
                       sgal_hist=sgal_opt_hist, sbins=sbins_opt,
                       model_title=model_title)
        logger.info('%s: %s OPT plotted' % (target, model))
        add_lines_LF(fig, axs, opt_gal, [opt_gal.trgb+0.5, opt_gal.trgb+2],
                     [opt_gal.model_plt_color, opt_gal.data_plt_color],
                     [ndata_agb_opt, nsim_agb_opt, ndata_rgb, nsim_rgb])  
        fig_name_kw['filter1'] = opt_gal.filter1
        fig_name_kw['filter2'] = opt_gal.filter2
        fig.savefig(figname % fig_name_kw, dpi=300,
                    bbox_to_inches='tight')
        fig.close()
    logger.info('IR %s %s %.3f $%.3f %.4f$ %.4f' % (opt_gal.target, model_short, nar_ratio_data, np.mean(nmodel), np.std(np.mean(nmodel)), mc_pct_diff))
    logger.info('OPT %s %s %.3f $%.3f %.4f$ %.4f' % (opt_gal.target, model_short, nar_ratio_data_opt, np.mean(nmodel_opt), np.std(np.mean(nmodel_opt)), mc_pct_diff_opt))

    print 'IR %s %s %.3f $%.3f\pm%.3f$ %.4f' % (opt_gal.target, model_short, nar_ratio_data, np.mean(nmodel), np.std(np.mean(nmodel)), mc_pct_diff)
    print 'OPT %s %s %.3f $%.3f\pm%.3f$ %.4f' % (opt_gal.target, model_short, nar_ratio_data_opt, np.mean(nmodel_opt), np.std(np.mean(nmodel_opt)), mc_pct_diff_opt)
    return


def multi_opt_rgb_nir_agb_ratio(filt1=None, filt2=None, targets=None, leo_ast=True, 
                          offsets=(1.5, 0.), models=None, maglims=None,
                          leo_method=False, leo_norm=False, make_plot=True,
                          xlim=None, ylim=None, xlim2=None, add_boxes=True,
                          color_hist=False,  plot_tpagb=False, **kwargs):

    ast_file_loc = snap_src + 'output/'

    # fix incoming args
    targets = load_targets(targets)
    if type(models) == str:
        models = [models]

    if make_plot is True:
        xlim, ylim, xlim2 = nir_cmd_plot_limits(xlim=xlim,
                                                ylim=ylim,
                                                xlim2=xlim2)

        plot_LF_kw = {'ylim': ylim, 'xlim': xlim, 'xlim2': xlim2,
                      'color_hist': color_hist, 'title': True}

    result_dict = {}
    print 'target model name data sim pct_dif'
    for target in targets:
        # used for agb_verts, will use box if not leo method, otherwise
        # will use all color space.
        if target in targets_leo_norm_ok():
            leo_method = True
        else:
            leo_method = False

        logger.info('working on %s' % target)

        # load data
        ir_gal = load_galaxy(target, band='ir')
        opt_gal = load_galaxy(target, band='opt')
        data_pts = np.column_stack((opt_gal.color, opt_gal.mag2))
        ir_data_pts = np.column_stack((ir_gal.color, ir_gal.mag2))

        # make the data histograms
        gal_opt_hist, gal_opt_bins = hist_it_up(opt_gal.mag2)
        gal_ir_hist, gal_ir_bins = hist_it_up(ir_gal.mag2)

        # get the rgb and agb polygons
        rgb_poly = get_opt_rgb_polygons(target)

        agb_verts = load_agb_verts(ir_gal, leo_method=leo_method)
        opt_agb_verts = load_agb_verts(opt_gal, leo_method=leo_method)

        if type(rgb_poly) == int:
            magdim = opt_gal.trgb + 2
            magbright = opt_gal.trgb + .5
            colmin = np.min(opt_gal.color)
            colmax = np.max(opt_gal.color)
            rgb_verts = np.array([[colmin, magdim],
                              [colmin, magbright],
                              [colmax, magbright],
                              [colmax, magdim],
                              [colmin, magdim]])       
        else:
            rgb_Color = rgb_poly[:, 0]
            rgb_Mag2 = rgb_poly[:, 1]
            rgb_Mag1 = rgb_Color + rgb_Mag2
            Mag2mag_kw = {'Av': opt_gal.Av, 'dmod': opt_gal.dmod}
            rmag1 = rsp.astronomy_utils.Mag2mag(rgb_Mag1, opt_gal.filter1, opt_gal.photsys,
                                                **Mag2mag_kw)

            rmag2 = rsp.astronomy_utils.Mag2mag(rgb_Mag2, opt_gal.filter2, opt_gal.photsys,
                                                **Mag2mag_kw)
            rgb_color = rmag1 - rmag2
            rgb_verts = np.column_stack([rgb_color, rmag2])

        # rgb and agb inds in the data
        rgb_in_data, = np.nonzero(nxutils.points_inside_poly(data_pts, rgb_verts))
        agb_in_data, = np.nonzero(nxutils.points_inside_poly(ir_data_pts, agb_verts))
        agb_in_data_opt, = np.nonzero(nxutils.points_inside_poly(data_pts, opt_agb_verts))

        # the nagb/nrgb ratio in the data
        ndata_rgb = float(len(rgb_in_data))                 
        ndata_agb = float(len(agb_in_data))
        nar_ratio_data = ndata_agb / ndata_rgb
        ndata_agb_opt = float(len(agb_in_data_opt))
        nar_ratio_data_opt = ndata_agb_opt / ndata_rgb

        # find the four filter ast file
        find_file = {'filter1': ir_gal.filter1, 
                     'filter2': ir_gal.filter2,
                     'filter3': opt_gal.filter1,
                     'filter4': opt_gal.filter2,
                     'target': opt_gal.target}
        figname = '%(target)s_%(model_name)s_%(filter1)s_%(filter2)s_LF.pdf'

        pool = multiprocessing.Pool()
        res = []
        for model in models:
            res.append(pool.apply_async(model_ratio_once, 
                                       (target, model, ast_file_loc, rgb_verts,
                                        ndata_rgb, agb_verts, opt_agb_verts,
                                        nar_ratio_data, nar_ratio_data_opt,
                                        ir_gal, opt_gal, find_file, figname, 
                                        ndata_agb, plot_LF_kw, ir_data_pts,
                                        gal_ir_hist, gal_ir_bins, gal_opt_hist,
                                        gal_opt_bins, data_pts, ndata_agb_opt),))

        for r in res:
            r.get()
    return


def count_uncert_ratio(numerator, denominator):
    '''
    combining poisson error for taking a ratio
    '''
    n = float(numerator)
    d = float(denominator)
    return (n / d) * (1./np.sqrt(n) + 1./np.sqrt(d))
 

def opt_rgb_nir_agb_ratio(filt1=None, filt2=None, targets=None, leo_ast=True, 
                          offsets=(1.5, 0.), models=None, maglims=None,
                          leo_method=False, leo_norm=False, make_plot=True,
                          xlim=None, ylim=None, xlim2=None, add_boxes=True,
                          color_hist=False,  plot_tpagb=False, norm_by_ir=False,
                          search_str=None, **kwargs):
    '''
    this will do either rgb norm in nir or opt...
    it's the right one to use, also hase uncertanties. it doesn't return 
    anything, but prints error to the table fmt (ish) and makes figures.    
    '''
    ast_file_loc = snap_src + 'output/'

    # fix incoming args
    targets = load_targets(targets)
    if type(models) == str:
        models = [models]

    if make_plot is True:
        xlim, ylim, xlim2 = nir_cmd_plot_limits(xlim=xlim,
                                                ylim=ylim,
                                                xlim2=xlim2)

        plot_LF_kw = {'ylim': ylim, 'xlim': xlim, 'xlim2': xlim2,
                      'color_hist': color_hist, 'title': True}

    print 'target model name data sim pct_dif'
    for target in targets:
        # used for agb_verts, will use box if not leo method, otherwise
        # will use all color space.
        if target in targets_leo_norm_ok():
            leo_method = True
        else:
            leo_method = False

        logger.info('working on %s' % target)

        # load data
        ir_gal = load_galaxy(target, band='ir')
        opt_gal = load_galaxy(target, band='opt')
        data_pts = np.column_stack((opt_gal.color, opt_gal.mag2))
        ir_data_pts = np.column_stack((ir_gal.color, ir_gal.mag2))

        # make the data histograms
        gal_opt_hist, gal_opt_bins = hist_it_up(opt_gal.mag2)
        gal_ir_hist, gal_ir_bins = hist_it_up(ir_gal.mag2)

        # get the rgb and agb polygons
        ir_agb_verts = load_agb_verts(ir_gal, leo_method=leo_method)
        opt_agb_verts = load_agb_verts(opt_gal, leo_method=leo_method)

        ir_gal.maglims = [ir_gal.trgb, ir_gal.trgb+1.5]
        rgb_verts_ir = get_ir_rgb_polygons(ir_gal, leo_method=leo_method)
        rgb_poly_opt = get_opt_rgb_polygons(target)

        if type(rgb_poly_opt) == int:
            magdim = opt_gal.trgb + 2
            magbright = opt_gal.trgb + .5
            colmin = np.min(opt_gal.color)
            colmax = np.max(opt_gal.color)
            rgb_verts_opt = np.array([[colmin, magdim],
                              [colmin, magbright],
                              [colmax, magbright],
                              [colmax, magdim],
                              [colmin, magdim]])       
        else:
            rgb_Color = rgb_poly_opt[:, 0]
            rgb_Mag2 = rgb_poly_opt[:, 1]
            rgb_Mag1 = rgb_Color + rgb_Mag2
            Mag2mag_kw = {'Av': opt_gal.Av, 'dmod': opt_gal.dmod}
            rmag1 = rsp.astronomy_utils.Mag2mag(rgb_Mag1, opt_gal.filter1, opt_gal.photsys,
                                                **Mag2mag_kw)

            rmag2 = rsp.astronomy_utils.Mag2mag(rgb_Mag2, opt_gal.filter2, opt_gal.photsys,
                                                **Mag2mag_kw)
            rgb_color = rmag1 - rmag2
            rgb_verts_opt = np.column_stack([rgb_color, rmag2])

        # rgb and agb inds in the data
        rgb_in_data_opt, = np.nonzero(nxutils.points_inside_poly(data_pts, rgb_verts_opt))
        rgb_in_data_ir, = np.nonzero(nxutils.points_inside_poly(ir_data_pts, rgb_verts_ir))
        agb_in_data_ir, = np.nonzero(nxutils.points_inside_poly(ir_data_pts, ir_agb_verts))
        agb_in_data_opt, = np.nonzero(nxutils.points_inside_poly(data_pts, opt_agb_verts))

        # the nagb/nrgb ratio in the data as well as errors
        ndata_rgb_opt = float(len(rgb_in_data_opt))                 
        ndata_rgb_ir = float(len(rgb_in_data_ir))
        
        ndata_agb_ir = float(len(agb_in_data_ir))
        if norm_by_ir is True:
            nar_ratio_data_ir = ndata_agb_ir / ndata_rgb_ir        
            nar_ratio_data_ir_err = count_uncert_ratio(ndata_agb_ir, ndata_rgb_ir)
        else:
            nar_ratio_data_ir = ndata_agb_ir / ndata_rgb_opt
            nar_ratio_data_ir_err = count_uncert_ratio(ndata_agb_ir, ndata_rgb_opt)
            
        ndata_agb_opt = float(len(agb_in_data_opt))
        nar_ratio_data_opt = ndata_agb_opt / ndata_rgb_opt

        nar_ratio_data_opt_err = count_uncert_ratio(ndata_agb_opt, ndata_rgb_opt)

        # find the four filter ast file
        find_file = {'filter1': ir_gal.filter1, 
                     'filter2': ir_gal.filter2,
                     'filter3': opt_gal.filter1,
                     'filter4': opt_gal.filter2,
                     'target': opt_gal.target}

        # printing fmts
        # will be:       band target model data, err model, err fdiff, ferr
        pfmt = '%s %s %s & $%.3f\pm%.3f$ & $%.3f\pm%.3f$ & $%.3f\pm%.3f$'
        if norm_by_ir is True:
            figname = '%(target)s_%(model_name)s_%(filter1)s_%(filter2)s_LF.pdf'
        else:
            figname = '%(target)s_%(model_name)s_%(filter1)s_%(filter2)s_opt_norm_LF.pdf'
        for model in models:
            logger.info('%s: %s' % (target, model))
            model_short = model.replace('cmd_input_', '').lower()

            # load four filter ast file
            find_file['model'] = model_short 
            if search_str is None:
                search_str = 'ast_%(filter1)s_%(filter2)s_%(filter3)s_%(filter4)s_output_%(target)s_%(model)s' % find_file
            sim_file, = rsp.fileIO.get_files(ast_file_loc, search_str)

            sgal = rsp.Galaxies.simgalaxy(sim_file, filter1='F110W', filter2='F160W')
            sgal.load_ast_corrections()

            # so i dunno why, but filter3 is 814, filter4 is the other guy.
            assert sgal.filter3 == 'F814W', 'uh uh.... %s' % sgal.filter3 
            
            opt_sim_pts = np.column_stack((sgal.ast_mag4 - sgal.ast_mag3,
                                           sgal.ast_mag3))
            ir_sim_pts = np.column_stack((sgal.ast_mag1 - sgal.ast_mag2,
                                          sgal.ast_mag2))
            
            # make normalization from the full simulation
            if norm_by_ir is True:
                rgb_in_sim, = np.nonzero(nxutils.points_inside_poly(ir_sim_pts,
                                                                    rgb_verts_ir))
                nsim_rgb_ir = float(len(rgb_in_sim))
                normalization = ndata_rgb_ir / nsim_rgb_ir
            else:
                rgb_in_sim, = np.nonzero(nxutils.points_inside_poly(opt_sim_pts,
                                                                    rgb_verts_opt))
                nsim_rgb_opt = float(len(rgb_in_sim))
                normalization = ndata_rgb_opt / nsim_rgb_opt

            # random sample the data distribution (all the same length)
            rands = np.random.random(len(sgal.ast_mag4))
            ind, = np.nonzero(rands < normalization)

            # scale the simulation
            norm_sim_pts_opt = opt_sim_pts[ind]
            norm_sim_pts_ir = ir_sim_pts[ind]

            # number of agb and rgb stars in the scaled simulation
            sim_agb_norm_ir, = np.nonzero(nxutils.points_inside_poly(norm_sim_pts_ir,
                                                                     ir_agb_verts))
            sim_rgb_norm_opt, = np.nonzero(nxutils.points_inside_poly(norm_sim_pts_opt,
                                                                      rgb_verts_opt))

            sim_rgb_norm_ir, = np.nonzero(nxutils.points_inside_poly(norm_sim_pts_ir,
                                                                     rgb_verts_ir))

            sim_agb_norm_opt, = np.nonzero(nxutils.points_inside_poly(norm_sim_pts_opt,
                                                                      opt_agb_verts))

            # for simulation's nagb/nrgb ratio
            nsim_agb_opt = float(len(sim_agb_norm_opt))
            nsim_agb_ir = float(len(sim_agb_norm_ir))
            nsim_rgb_ir = float(len(sim_rgb_norm_ir))
            nsim_rgb_opt = float(len(sim_rgb_norm_opt))
            
            # if we are norm by optical use that instead of ir...
            if norm_by_ir is True:
                nsim_rgb_ir = nsim_rgb_opt

            # nagb/nrgb in the scaled simulation
            nar_ratio_sim_ir = nsim_agb_ir / nsim_rgb_ir
            nar_ratio_sim_opt = nsim_agb_opt / nsim_rgb_opt

            # calculate errors
            nar_ratio_sim_ir_err = count_uncert_ratio(nsim_agb_ir, nsim_rgb_ir)
            nar_ratio_sim_opt_err = count_uncert_ratio(nsim_agb_opt, nsim_rgb_opt)
            
            # calculate fractional differences
            pct_diff_ir = (nar_ratio_sim_ir - nar_ratio_data_ir) / nar_ratio_data_ir
            pct_diff_opt = (nar_ratio_sim_opt - nar_ratio_data_opt) / nar_ratio_data_opt

            # propagate uncertainties
            pct_diff_ir_err = pct_diff_ir * (nar_ratio_sim_ir_err/nar_ratio_sim_ir + nar_ratio_data_ir_err/nar_ratio_data_ir)
            pct_diff_opt_err = pct_diff_opt * (nar_ratio_sim_opt_err/nar_ratio_sim_opt + nar_ratio_data_opt_err/nar_ratio_data_opt)
            
            # make histogram, bins for plotting.
            sgal_ir_hist, sbins_ir = hist_it_up(norm_sim_pts_ir[:, 1])
            sgal_opt_hist, sbins_opt = hist_it_up(norm_sim_pts_opt[:, 1])

            # plot the IR LF
            model_title = translate_model_name(model)
            plot_LF_kw['model_title'] = model_title
            fig, axs, _ = ir_gal.plot_LF(ir_data_pts[:,0], ir_data_pts[:,1],
                         norm_sim_pts_ir[:, 0], norm_sim_pts_ir[:, 1],
                         ir_gal.filter1, ir_gal.filter2,
                         gal_hist=gal_ir_hist, bins=gal_ir_bins,
                         sgal_hist=sgal_ir_hist, sbins=sbins_ir,
                         **plot_LF_kw)
            # add the numbers of each population to the plots
            add_lines_LF(fig, axs, ir_gal, [ir_gal.trgb, ir_gal.trgb+2],
                         [ir_gal.model_plt_color, ir_gal.data_plt_color],
                         [ndata_agb_ir, nsim_agb_ir, ndata_rgb_ir, nsim_rgb_ir])  
            
            # save fig
            fig_name_kw = {'target': opt_gal.target,
                           'model_name': model.replace('cmd_input_', '').replace('.dat', ''),
                           'filter1': ir_gal.filter1,
                           'filter2': ir_gal.filter2}
            plt.savefig(figname % fig_name_kw, dpi=300,
                        bbox_to_inches='tight')
            
            # print table to screen
            print pfmt % ('IR', ir_gal.target, model_short, nar_ratio_data_ir,
                          nar_ratio_data_ir_err, nar_ratio_sim_ir, nar_ratio_sim_ir_err,
                          pct_diff_ir, pct_diff_ir_err)
            plot_comp_LF()
            # repeat above output steps in OPT if this is a G10 galaxy
            if target in gi10_overlap():
                # if it's part of G10 plot the OPT LF
                fig, axs, _ = opt_gal.plot_LF(data_pts[:,0], data_pts[:,1],
                               norm_sim_pts_opt[:, 0], norm_sim_pts_opt[:, 1],
                               opt_gal.filter1, opt_gal.filter2,
                               gal_hist=gal_opt_hist, bins=gal_opt_bins,
                               sgal_hist=sgal_opt_hist, sbins=sbins_opt,
                               model_title=model_title)
                add_lines_LF(fig, axs, opt_gal, [opt_gal.trgb+0.5, opt_gal.trgb+2],
                             [opt_gal.model_plt_color, opt_gal.data_plt_color],
                             [ndata_agb_opt, nsim_agb_opt, ndata_rgb_opt, nsim_rgb_opt])  
                fig_name_kw['filter1'] = opt_gal.filter1
                fig_name_kw['filter2'] = opt_gal.filter2
                plt.savefig(figname % fig_name_kw, dpi=300,
                            bbox_to_inches='tight')
                print pfmt % ('OPT', opt_gal.target, model_short,
                              nar_ratio_data_opt, nar_ratio_data_opt_err,
                              nar_ratio_sim_opt, nar_ratio_sim_opt_err,
                              pct_diff_opt, pct_diff_opt_err)
                              


    return
    

def plot_comp_LF(models, targets, norm_by_ir=False ):
    # fix incoming args
    targets = load_targets(targets)
    if type(models) == str:
        models = [models]

    young_color = 'navy'
    old_color = 'darkred'
    for model in models:
        model_short = model.replace('cmd_input_', '').lower()
        fig, (ax1, ax2, ax3) = plt.subplots(nrows=3, figsize=(10,12))
        resids = []
        pct_diffs = []
        model_name = translate_model_name(model)
        for target in targets:
            leo_method = target in targets_leo_norm_ok()

            # load data
            ir_gal = load_galaxy(target, band='ir')
            opt_gal = load_galaxy(target, band='opt')
            data_pts = np.column_stack((opt_gal.color, opt_gal.mag2))
            ir_data_pts = np.column_stack((ir_gal.color, ir_gal.mag2))

            ir_gal.maglims = [ir_gal.trgb, ir_gal.trgb+1.5]

            rgb_verts_ir = get_ir_rgb_polygons(ir_gal, leo_method=leo_method)
            rgb_verts_opt = get_rgb_verts_opt(opt_gal)

            # get the rgb and agb polygons
            ir_agb_verts = load_agb_verts(ir_gal, leo_method=leo_method)
            opt_agb_verts = load_agb_verts(opt_gal, leo_method=leo_method)


            # rgb and agb inds in the data
            rgb_in_data_opt, = np.nonzero(nxutils.points_inside_poly(data_pts, rgb_verts_opt))
            rgb_in_data_ir, = np.nonzero(nxutils.points_inside_poly(ir_data_pts, rgb_verts_ir))
            agb_in_data_ir, = np.nonzero(nxutils.points_inside_poly(ir_data_pts, ir_agb_verts))
            agb_in_data_opt, = np.nonzero(nxutils.points_inside_poly(data_pts, opt_agb_verts))

            # the nagb/nrgb ratio in the data as well as errors
            ndata_rgb_opt = float(len(rgb_in_data_opt))                 
            ndata_rgb_ir = float(len(rgb_in_data_ir))
        
            ndata_agb_ir = float(len(agb_in_data_ir))
            if norm_by_ir is True:
                nar_ratio_data_ir = ndata_agb_ir / ndata_rgb_ir        
                nar_ratio_data_ir_err = count_uncert_ratio(ndata_agb_ir, ndata_rgb_ir)
            else:
                nar_ratio_data_ir = ndata_agb_ir / ndata_rgb_opt
                nar_ratio_data_ir_err = count_uncert_ratio(ndata_agb_ir, ndata_rgb_opt)
            
            ndata_agb_opt = float(len(agb_in_data_opt))
            nar_ratio_data_opt = ndata_agb_opt / ndata_rgb_opt

            nar_ratio_data_opt_err = count_uncert_ratio(ndata_agb_opt, ndata_rgb_opt)

            # find the four filter ast file
            find_file = {'filter1': ir_gal.filter1, 
                         'filter2': ir_gal.filter2,
                         'filter3': opt_gal.filter1,
                         'filter4': opt_gal.filter2,
                         'target': opt_gal.target}

            # printing fmts
            # will be:       band target model data, err model, err fdiff, ferr
            pfmt = '%s %s %s & $%.3f\pm%.3f$ & $%.3f\pm%.3f$ & $%.3f\pm%.3f$'
            if norm_by_ir is True:
                figname = '%(target)s_%(model_name)s_%(filter1)s_%(filter2)s_LF.pdf'
            else:
                figname = '%(target)s_%(model_name)s_%(filter1)s_%(filter2)s_opt_norm_LF.pdf'
            sim_file = load_ast_file(ir_gal, model)
            sgal = rsp.Galaxies.simgalaxy(sim_file, filter1='F110W', filter2='F160W')
            sgal.load_ast_corrections()

            ir_sim_pts = np.column_stack((sgal.ast_mag1 - sgal.ast_mag2,
                                          sgal.ast_mag2))
        
            # make normalization from the full simulation
            rgb_in_sim, = np.nonzero(nxutils.points_inside_poly(ir_sim_pts,
                                                                rgb_verts_ir))
            nsim_rgb_ir = float(len(rgb_in_sim))
            normalization = ndata_rgb_ir / nsim_rgb_ir

            nsim_rgb_opt = float(len(rgb_in_sim))
            normalization = ndata_rgb_opt / nsim_rgb_opt

            # random sample the data distribution (all the same length)
            rands = np.random.random(len(sgal.ast_mag2))
            ind, = np.nonzero(rands < normalization)
    
            # scale the simulation
            norm_sim_pts_ir = ir_sim_pts[ind]

            # number of agb and rgb stars in the scaled simulation
            sim_agb_norm_ir, = np.nonzero(nxutils.points_inside_poly(norm_sim_pts_ir,
                                                                     ir_agb_verts))

            sim_rgb_norm_ir, = np.nonzero(nxutils.points_inside_poly(norm_sim_pts_ir,
                                                                     rgb_verts_ir))


            # for simulation's nagb/nrgb ratio
            nsim_agb_ir = float(len(sim_agb_norm_ir))
            nsim_rgb_ir = float(len(sim_rgb_norm_ir))
        
            # if we are norm by optical use that instead of ir...

            # nagb/nrgb in the scaled simulation
            nar_ratio_sim_ir = nsim_agb_ir / nsim_rgb_ir

            # calculate errors
            nar_ratio_sim_ir_err = count_uncert_ratio(nsim_agb_ir, nsim_rgb_ir)
        
            # calculate fractional differences
            pct_diff_ir = (nar_ratio_sim_ir - nar_ratio_data_ir) / nar_ratio_data_ir

            # propagate uncertainties
            pct_diff_ir_err = pct_diff_ir * (nar_ratio_sim_ir_err/nar_ratio_sim_ir + nar_ratio_data_ir_err/nar_ratio_data_ir)
        
            # make histogram, bins for plotting.
            sgal_ir_hist, sbins_ir = hist_it_up(norm_sim_pts_ir[:, 1])

            # make the data histograms --- use the bins from above!
            gal_ir_hist = np.histogram(ir_gal.mag2, bins=sbins_ir)[0]

            bins = rsp.astronomy_utils.mag2Mag(sbins_ir[1:], 'F160W',
                                               'wfc3snap',
                                               **{'target': ir_gal.target,
                                                  'filter1': ir_gal.filter1})
                                                  
            # set to get all mags... should be Trgb if you want to slice.
            comp_inds = np.nonzero(bins < ir_gal.Trgb)

            pct_diff = (sgal_ir_hist - gal_ir_hist) / gal_ir_hist
            resid = sgal_ir_hist - gal_ir_hist
            resids.append(np.mean(resid[comp_inds]))
            pct_diffs.append(np.mean(pct_diff[comp_inds]))
            # 3 panel
            # top: # vs magnitude
            # middle: residual model-data  vs mag
            # bottom: %diff vs mag
            if target in ['DDO71', 'DDO78', 'DDO82', 'SCL-DE1']:
                col = old_color
            else:
                col = young_color
    
            label = '$%s$' % ir_gal.target.replace('-', '\!-\!')
            plt_kw = {'ls': 'steps', 'color': col, 'lw': 2}
            #[ax.vlines(smg.gal.Trgb, *ax.get_ylim(), color=cols[j]) for ax in [ax1, ax2, ax3]]
            ax1.semilogy(bins[comp_inds], sgal_ir_hist[comp_inds], label=label, **plt_kw)
            ax2.plot(bins[comp_inds], pct_diff[comp_inds], **plt_kw)
            ax3.plot(bins[comp_inds], resid[comp_inds], **plt_kw)
            #plt_kw['ls'] += '--'
            #plt_kw['lw'] = 1
            #ax1.semilogy(bins[comp_inds], gal_ir_hist[comp_inds], **plt_kw)
            print ir_gal.target, np.mean(resid[comp_inds]), np.std(resid[comp_inds]), np.mean(pct_diff[comp_inds]), np.std(pct_diff[comp_inds])
        print model_name, np.mean(resids), np.std(resids), np.mean(pct_diffs), np.std(pct_diffs)
        ax1.text(0.75, 0.85, '$%s$' % model_name, fontsize=20, transform=ax1.transAxes)
        ax2.text(0.75, 0.85, '${\\rm Frac.\ Difference}$', fontsize=20,
                 transform=ax2.transAxes)        
        ax3.text(0.75, 0.85, '${\\rm Residual}$', fontsize=20,
                 transform=ax3.transAxes)

        #ax1.legend(loc=0, numpoints=1, frameon=False)
        ax1.set_ylabel('$\#/{\\rm mag}$', fontsize=20)
        ax2.set_ylabel('$(N_{\\rm model}-N_{\\rm data})/N_{\\rm data}$', fontsize=20)
        ax3.set_ylabel('$N_{\\rm model}-N_{\\rm data}$', fontsize=20)
        [ax.set_xlabel('$F160W$', fontsize=20) for ax in [ax1, ax2, ax3]]
        plt.tick_params(labelsize=16)
        ax1.set_ylim(10, ax1.get_ylim()[1])
        ax2.set_ylim(-2, 50)
        ax3.set_ylim(-100, 600)
        [ax.set_xlim(-5.5, -9) for ax in [ax1, ax2, ax3]]
        plt.subplots_adjust(hspace=.25, top=.95)
        mname = model.replace('.dat','').split('_')[-1]
        plt.savefig('comp_lfs_%s.png' % mname, dpi=300)

    return


def get_rgb_verts_opt(opt_gal):
    rgb_poly_opt = get_opt_rgb_polygons(opt_gal.target)

    if type(rgb_poly_opt) == int:
        magdim = opt_gal.trgb + 2
        magbright = opt_gal.trgb + .5
        colmin = np.min(opt_gal.color)
        colmax = np.max(opt_gal.color)
        rgb_verts_opt = np.array([[colmin, magdim],
                          [colmin, magbright],
                          [colmax, magbright],
                          [colmax, magdim],
                          [colmin, magdim]])       
    else:
        rgb_Color = rgb_poly_opt[:, 0]
        rgb_Mag2 = rgb_poly_opt[:, 1]
        rgb_Mag1 = rgb_Color + rgb_Mag2
        Mag2mag_kw = {'Av': opt_gal.Av, 'dmod': opt_gal.dmod}
        rmag1 = rsp.astronomy_utils.Mag2mag(rgb_Mag1, opt_gal.filter1, opt_gal.photsys,
                                            **Mag2mag_kw)

        rmag2 = rsp.astronomy_utils.Mag2mag(rgb_Mag2, opt_gal.filter2, opt_gal.photsys,
                                            **Mag2mag_kw)
        rgb_color = rmag1 - rmag2
        rgb_verts_opt = np.column_stack([rgb_color, rmag2])
    return rgb_verts_opt


def add_lines_LF(fig, axs, gal, maglims, colors, vals):
    '''
    must have attributes sgal, gal nbrighter
    '''
    j = 0
    line_on_it_kw = {'annotate': 0, 'ls': '-'}
    for i, maglim in enumerate(maglims):
        # lines and numbers on plots
        for ax, col in zip(axs[:2], colors):
            gal.put_a_line_on_it(ax, maglim, color=col, **line_on_it_kw)
            gal.annotate_cmd(ax, maglim, '$%i$' % vals[j], text_kw={'color': col})
            j += 1
    return fig, axs


def hist_it_up(mag2, res=0.1):
    # do fine binning and hist
    bins = np.arange(np.nanmin(mag2), np.nanmax(mag2), res)
    hist = np.histogram(mag2, bins=bins)[0]
    # drop the too fine bins
    binds = spread_bins(hist)
    # return the hist, bins.
    return np.histogram(mag2, bins=bins[binds])


def spread_bins(hist, threash=10):
    '''
    goes through in index order of hist and returns indices such that
    each bin will add to at least threash.
    Returns the indices of the new bins to use (should then re-hist)
    ex:
    bins is something set that works well for large densities
    h = np.histogram(mag2, bins=bins)[0]
    binds = spread_bins(h)
    hist, bins = np.histogram(mag2, bins=bins[binds])
    '''
    b = []
    i = 1
    j = 0
    while i < len(hist):
        while np.sum(hist[j:i]) < threash:
            i += 1
            #print j, i, len(hist)
            if i > len(hist):
                break
        j = i
        b.append(j)
    return np.array(b[:-1])


def ir_rgb_agb_ratio(renormalize=False, filt1=None, filt2=None, band=None,
                     targets=None, run_trilegal=False, leo_ast=True, 
                     offsets=(1.5, 0.), models=None, maglims=None,
                     leo_method=False, leo_norm=False, make_plot=True,
                     xlim=None, ylim=None, xlim2=None, add_boxes=True,
                     color_hist=False,  plot_tpagb=False, use_opt_rgb=False,
                     sfr_dir='default', **kwargs):
    smgs = []
    targets = load_targets(targets)

    if type(models) == str:
        models = [models]

    norm_sim_kw = {'offsets': offsets, 'band': band, 'leo_ast': leo_ast,
                   'run_trilegal': run_trilegal, 'leo_method': leo_method,
                   'leo_norm': leo_norm, 'use_opt_rgb': use_opt_rgb,
                   'sfr_dir': sfr_dir}

    if make_plot is True:
        if band == 'ir':
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
    models_loc = snap_src + 'models/match/'
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
        model_name = '$\dot{M}_{G10_b}$'
    elif 'mar13' in model.lower():
        model_name = '$\dot{M}_{G10}$'
    elif 'oct13' in model.lower():
        model_name = '$\dot{M}_{M13}$'
    return model_name


def fuck_you_match(gal, sgal, verts=False, inverse_verts=False, make_plot=False,
                   figname=None, dmag=0.1, dcol=0.05, modelB=None,
                   rgb_only=False, labels=None, skip_agb=False, band='ir'):

    if rgb_only is True:
        sgal.all_stages('RGB')
        inds = np.intersect1d(sgal.irgb, sgal.norm_inds)
    elif skip_agb is True:
        sgal.all_stages('TPAGB')
        inds = [i for i in sgal.norm_inds if i not in sgal.itpagb]
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
        elif skip_agb is True:
            gal.all_stages('TPAGB')
            ginds = [i for i in gal.norm_inds if i not in gal.itpagb]
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
            figname = 'pgpro_%s_%s_%s_%s%s.pdf' % (gal.target, model, gal.filter1, gal.filter2, extra)        

        model_name = translate_model_name(sgal.model)

        ZS = [data_hess[2], sim_hess[2], pct_dif, sig]
        xlim, ylim, _ = nir_cmd_plot_limits()
        if band == 'opt':
            xlim = (-0.5, 2.5)
            ylim = (28, 18)

        extent = [data_hess[0][0], data_hess[0][-1], data_hess[1][-1], data_hess[1][0]]
        if labels is None:
            labels = ['${\\rm %s}$' % gal.target, model_name, '${\\rm Frac.\ Difference}$', '$\chi^2=%.2f$' % chi2]
        if labels[2] is None:
            labels[2] = '${\\rm Frac.\ Difference}$'
        if labels[3] is None:
            labels[3] = '$\chi^2=%.2f$' % chi2
            
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
        print gal.target, gal.filter1, gal.filter2, Nrgb, len(inds)
    return Nrgbs


def run_fuck_you_match(targets=None, models=None, band=None, inputs=None,
                       modelB=None):
    '''
    does the chi2 test for full field, agb region, and not-agb region.
    to compare model to model, use modelB.
    
    this is called fuck you match because I realized at 3 am that 
    match smooths the bg cmd when doing stats tests, and was ruining 
    my results, so I had to write my own (which didn't take long, but
    still)
    '''
    inputs = inputs or {}
    if targets is None:
        targets = inputs['targets']
    if band is None:
        band = inputs['band']
    if models is None:
        models = inputs['models']

    targets = load_targets(targets)

    if type(models) == str:
        models = [models]

    figname_fmt = 'pgpro_%s_%s_%s.pdf'

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
                model_name = translate_model_name(model)
                modelB_name = translate_model_name(modelB)
                if 'MAR13' in model:
                    model_name = '${\\rm PARSEC}$'
                labels = ['${\\rm Padova}$', model_name, None, None]

                #gal.target = modelB.replace('.dat', '').split('_')[-1]
            # entire cmd fit
            extra = 'full'
            print extra
            

            figname = figname_fmt % (target, model.replace('.dat',''), extra)
            full_chi2 = fuck_you_match(gal, sgal, figname=figname, modelB=modelB,
                                       make_plot=True, labels=labels, band=band)
            
            if band == 'opt':
                print '%s %s %.2f' % (target, model.replace('.dat', '').split('_')[-1], full_chi2)
                continue
            # include agb only
            extra = 'agb'
            print extra
            figname = figname_fmt % (target, model.replace('.dat',''), extra)
            agb_chi2 = fuck_you_match(gal, sgal, figname=figname, verts=True,
                                      make_plot=True, modelB=modelB, labels=labels)
            # exclude agb
            extra = 'not_agb'
            print extra
            figname = figname_fmt % (target, model.replace('.dat',''), extra)
            part_chi2 = fuck_you_match(gal, sgal, figname=figname, modelB=modelB,
                                       inverse_verts=True, make_plot=True, labels=labels)

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
        bg_loc = snap_src + 'models/match/bg/'
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
        mname = models[i].replace('.dat','').split('_')[-1]
        plt.savefig('comp_lfs_%s.png' % mname, dpi=300)

    do_a_plot_that_wastes_time = True
    if do_a_plot_that_wastes_time is True:
        for i in range(len(models)):
            for ycol in ['age', 'mass']:
                model_name = translate_model_name(models[i])
                fig, ax  = plt.subplots(figsize=(10, 10))

                divider = make_axes_locatable(ax)
                axt = divider.append_axes("top", 1.5, pad=0.1, sharex=ax)
                axr = divider.append_axes("right", 1.5, pad=0.1, sharey=ax)

                ax.text(1.12, 1.12, model_name, fontsize=20,
                        transform=ax.transAxes)
                for j, smg in enumerate(np.array(smgs)[np.where(imod==i)]):
                    plt_kw = {'histtype': 'step', 'color': cols[j], 'lw': 2}
                    smg.sgal.all_stages('TPAGB')
                
                    points = np.column_stack((smg.sgal.ast_color,
                                              smg.sgal.ast_mag2))
                    inregion, = np.nonzero(nxutils.points_inside_poly(points,
                                                                      smg.agb_verts))

                    tpagb_inds = list(set(smg.sgal.norm_inds) &
                                      set(smg.sgal.itpagb) &
                                      set(inregion))

                    poop = np.intersect1d(smg.sgal.norm_inds, inregion) 
                    not_tpagb = np.array([k for k in poop
                                          if k not in smg.sgal.itpagb])
                    label='$%s$' % smg.gal.target.replace('-','\!-\!')

                    m2M_kw = {'target': smg.gal.target,
                              'filter1': smg.gal.filter1}
                    mag2 = rsp.astronomy_utils.mag2Mag(smg.sgal.ast_mag2,
                                                       'F160W', 'wfc3snap',
                                                        **m2M_kw)
                    x = mag2[tpagb_inds]
                    xx = mag2[not_tpagb]
                    if ycol == 'age':
                        ydata = smg.sgal.data.get_col('logAge')
                    if ycol == 'mass':
                        ydata = smg.sgal.data.get_col('m_ini')

                    y = ydata[tpagb_inds]
                    yy = ydata[not_tpagb]
    
                    ax.plot(x, y, 'o', ms=5, alpha=.5, color=cols[j],
                            mec=cols[j], label=label)
            
                    ax.set_xlim(-10., -5.5)
                    if ycol == 'mass':
                        ax.set_ylim(0.5, 6.)
                    if ycol == 'age':
                        ax.set_ylim(10.5, 7.5)
                    # make some labels invisible
                    plt.setp(axt.get_xticklabels() + axr.get_yticklabels(),
                             visible=False)

                    bins = int(np.round(np.sqrt(len(x))))
                    axt.hist(x, bins=bins, **plt_kw)
                    axr.hist(y, bins=bins, log=True, orientation='horizontal',
                             **plt_kw)
            
                ax.set_xlabel('$F160W$', fontsize=20)
                if ycol == 'mass':
                    ax.set_ylabel('$M (M_\odot)$', fontsize=20)
                if ycol == 'age':
                    ax.set_ylabel('${\log \\rm Age\ (yr)}$', fontsize=20)

                for tl in axt.get_xticklabels():
                    tl.set_visible(False)

                for tl in axr.get_yticklabels():
                    tl.set_visible(False)

                ax.legend(loc=0, numpoints=1, frameon=False)
                mname = models[i].replace('.dat','').split('_')[-1]
                plt.savefig('%s_mag_%s.png' % (ycol, mname), dpi=300)


def sfh_plots(sfh_loc=None, targets='paper1', make_plot=False,
                 ext='.dat', dan_fmt=False):
    if sfh_loc is None:
        sfh_loc = snap_src + 'data/sfh/'
        sfh_loc = '/Users/phil/research/Italy/WFC3SNAP/PHIL/SFRfiles/noAGB/fullreszctmps/'
        
    #bmap = brewer2mpl.get_map('Paired', 'Qualitative', 9)
    #cols = bmap.mpl_colors

    fig, ax = plt.subplots()
    for j, target in enumerate(targets):
        sfh_file, = rsp.fileIO.get_files(sfh_loc, '*%s*%s' % (target, ext))
        print sfh_file
        if not target in sfh_file:
            continue
        if dan_fmt is False:
            x = sfh_file.split('_')[2]
        else:
            x = os.path.split(sfh_file)[1].split('.')[0]
        try:
            a, s, z = np.loadtxt(sfh_file, unpack=True)
        except ValueError:
            to, tf, s, mh = np.genfromtxt(sfh_file, skip_header=6, skip_footer=2, usecols=[0,1,3,6], unpack=True)
            a = (to + tf) / 2
            from ResolvedStellarPops.convertz import convertz
            z = np.array([convertz(mh=i)[1] for i in mh])

        ax.plot(a[::-1], np.cumsum(s[::-1])/np.sum(s), lw=3, label=target,
                color=col)
    ax.set_xlim(10.13, 6.6)
    ax.set_ylim(0, 1)
    ax.hlines(0.6, *ax.get_xlim(), color='black', lw=4)
    ax.set_ylabel('${\\rm Cumulative\ Star\ Formation}$', fontsize=20)
    ax.set_xlabel('$\log {\\rm Time\ (yrs\ ago)}$', fontsize=20)
    plt.tick_params(labelsize=16)
    ax.legend(loc=0, frameon=False)
    plt.savefig('cumsum_sfr.pdf', dpi=300)
    

def match_metals(sfh_loc=None, targets=None, make_plot=False,
                 ext='.dat', dan_fmt=False):
    '''
    print the sfh weighted metallicity from the match input file. 
    '''
    if sfh_loc is None:
        sfh_loc = snap_src + 'data/sfh/'
    sfh_files = rsp.fileIO.get_files(sfh_loc, '*%s' % ext)
    if targets is None:
        targets = sfh_files
    for sfh_file in sfh_files:
        print sfh_file
        if dan_fmt is False:
            target = sfh_file.split('_')[2]
        else:
            target = os.path.split(sfh_file)[1].split('.')[0]
        try:
            a, s, z = np.loadtxt(sfh_file, unpack=True)
        except ValueError:
            to, tf, s, mh = np.genfromtxt(sfh_file, skip_header=6, skip_footer=2, usecols=[0,1,3,6], unpack=True)
            a = (to + tf) / 2
            from ResolvedStellarPops.convertz import convertz
            z = np.array([convertz(mh=i)[1] for i in mh])

        o, = np.nonzero(np.log10(a*1e9) < 8.5)
        
        print '%.4f %.4f %.4f %s ' % (np.min(z), np.sum(z*s)/np.sum(s), np.max(z), target)
        #print '%.4f %.4f %.4f %s ' % (np.min(z[o]), np.sum(z[o]*s[o])/np.sum(s[o]), np.max(z[o]), sfh_file.split('_')[2])
        if make_plot is True:
            fig, ax = plt.subplots()
            ax.plot(a, s)
            ax.set_title(target)
            fig.savefig(target.replace('.dat','_sfh.png'))


def agb_logl_age():
    for i in range(len(models)):
        model_name = models[i].replace('.dat', '').split('_')[-1]
        agb_track_loc = research_path + 'AGBTracks/CAF09/S_%s/S12_Z0.002_Y0.252/' % model_name
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
    import fileIO
    tauss = []
    btauss = []
    bmap = brewer2mpl.get_map('Set1', 'Qualitative', 3)
    cols = bmap.mpl_colors
    fig, ax = plt.subplots()
    fig2, ax2 = plt.subplots()
    for i in range(len(models)):
        model_name = models[i].replace('.dat', '').split('_')[-1]
        agb_track_loc = research_path + 'AGBTracks/CAF09/S_%s/S12_Z0.002_Y0.252/' % model_name
        if not os.path.isdir(agb_track_loc) is True:
            print model_name, 'no agb tracks found'
        model_name = translate_model_name(models[i])
        agb_track_names = [os.path.join(agb_track_loc, a)
                           for a in os.listdir(agb_track_loc)
                           if a.startswith('agb_')]
        tracks = [fileIO.get_numeric_data(agb_track)
                  for agb_track in agb_track_names]
        masses = np.array([t.mass for t in tracks])
        logls = np.array([t.get_col('L_star') for t in tracks])
        brights = np.array([np.nonzero(logl > 3.4)[0] for logl in logls])
        #m_cs = np.array([t.get_col('M_c')[0] for t in tracks])
        #ax2.plot(masses, m_cs, lw=2, color='black')

        taus = np.array([np.sum(t.data_array['dt']) for t in tracks])
        btaus = np.array([np.sum(t.data_array['dt'][b]) for t,b in zip(tracks, brights)])
        tauss.append(taus)
        btauss.append(btaus)
        ax.plot(masses, taus/1e6, lw=3, label=model_name, color=cols[i])
        ax2.plot(masses, btaus/1e6, lw=3, label=model_name, color=cols[i])

    ax.fill_between(masses, tauss[0]/1e6, tauss[2]/1e6, alpha=0.1, color='grey')
    ax2.fill_between(masses, btauss[0]/1e6, btauss[2]/1e6, alpha=0.1, color='grey')


    for ax in [ax, ax2]:
        ax.set_xlabel('${\\rm Initial\ Mass\ (M_\odot)}$', fontsize=20)
        ax.set_ylabel('${\\rm Lifetime\ (Myr)}$', fontsize=20)
        ax.legend(loc=0, frameon=False)
        ax.set_xlim(0, 5)
        ax.set_ylim(0, 7)
    #ax2.set_ylabel('${\\rm Pre\!-\!Flash\ Core\ Mass\ (M_\odot)}$', fontsize=20)
    fig.savefig('tpagb_lifetime.png', dpi=300)
    fig2.savefig('tpagb_lifetime_bright.png', dpi=300)


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


def fuckshitballs():

    targets = load_targets('gi10')
    opt_fakes = [galaxy_tests.get_fake_files(t, band='opt') for t in targets]

    filter2 = 'F814W'
    filter1s = ['F475W', 'F606W', 'F606W']

    sgals = [rsp.Galaxies.simgalaxy(tri_outs[i], filter1=filter1s[i], 
                                    filter2=filter2, photsys='acs_wfc')
             for i in range(len(targets))]

    tri_outs =[rsp.fileIO.get_files(snap_src + 'output/','output_*%s*_opt.dat' % t)[0] for t in targets]
    sgals = [rsp.Galaxies.simgalaxy(tri_outs[i], filter1=filter1s[i],
                                    filter2=filter2, photsys='acs_wfc')
             for i in range(len(targets))]


    outfiles = [os.path.join(sgals[i].base, 'ast_%s' % sgals[i].name)
                for i in range(len(sgals))]


    [rsp.Galaxies.ast_correct_trilegal_sim(sgals[i], fake_file=fakes[i],
                                           outfile=outfiles[i],
                                           leo_method=False)
     for i in range(len(sgals))]

    models = 'cmd_input_gi10_rev_old_tracks.dat'

    opt_fakes = [galaxy_tests.get_fake_files(t, band='opt')
                 for t in targets]

    opt_gals = [galaxy_tests.load_galaxy(t, band='opt') for t in targets]

    ir_gals = [galaxy_tests.load_galaxy(t, band='ir') for t in targets]

    fourfilts = ['_'.join(np.concatenate([outfiles[i].split('_')[:3],
                                          ir_gals[i].filters,
                                          opt_gals[i].filters,
                                          outfiles[i].split('_')[3:]]))
                 for i in range(len(outfiles))]

    sgal2s = [rsp.Galaxies.simgalaxy(outfiles[i], filter1=opt_gals[i].filter1,
                                     filter2=opt_gals[i].filter2)
              for i in range(len(outfiles))]
    [rsp.Galaxies.ast_correct_trilegal_sim(sgal2s[i], fake_file=opt_fakes[i],
                                           outfile=fourfilts[i],
                                           leo_method=False)
     for i in range(len(sgal2s))]
    sgals[0].plot_cmd(sgals[0].color, sgals[0].mag2, scatter_off=True)
    search_string = 'ast_output_%(target)s_%(filter1)s_%(filter2)s_%(filter3)s_%(filter4)s_*opt.dat'

    galaxy_tests.opt_rgb_nir_agb_ratio(targets='gi10', models=models, search_string=search_string)


if __name__ == "__main__":
    #snap_src = research_path + 'SNAP'
    #models =  ['gi10_rev_old_tracks.dat', 'cmd_input_CAF09_S_APR13.dat',
    #           'cmd_input_CAF09_S_APR13VW93.dat', 'cmd_input_CAF09_S_MAR13.dat']
    #opt_rgb_nir_agb_ratio(targets='paper1', models=models, norm_by_ir=True)
    main(sys.argv[1])
else:
    global snap_src
    if 'Linux' in os.uname():
        snap_src = '/home/phil/research/TP-AGBcalib/SNAP'
    else:
        snap_src = research_path + 'SNAP'
