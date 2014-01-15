import ResolvedStellarPops as rsp
import numpy as np
import galaxy_tests
import sfh_tests
import os
import fileIO
import matplotlib.pyplot as plt
from TPAGBparams import research_path, snap_src
from astroML.stats import binned_statistic_2d
import matplotlib.gridspec as gridspec
import brewer2mpl
from matplotlib.colors import LogNorm

'''
This is a work in progress, which sorts of intro figures should I have for the paper??

'''


def plot_opt_hess(targets=None, fits_src='default', filter1='F606W'):
    band = 'opt'
    ylim = (0.5, -7)
    if filter1 == 'F110W':
        band = 'ir'
        fits_src = 'default'
        ylim = (-2, -8)
    targets = galaxy_tests.load_targets(targets)
    galss = rsp.Galaxies.galaxies([galaxy_tests.load_galaxy(t, band=band,
                                                           fits_src=fits_src)
                                  for t in targets])

    gals = rsp.Galaxies.galaxies(galss.select_on_key('filter1', filter1))

    gals.squish('Color', 'Mag2', 'Trgb')

    mean_trgb = np.mean(gals.Trgbs)
    offset = gals.Trgbs - np.mean(gals.Trgbs)
    gals.Mag2o = np.concatenate([g.Mag2 - offset[i] for i, g in enumerate(gals.galaxies)])
    gals.Colorso = gals.Colors + gals.Mag2s - gals.Mag2o
    galshess =  rsp.astronomy_utils.hess(gals.Colorso, gals.Mag2o, 0.1, cbinsize=0.05)

    #N, xedges, yedges = binned_statistic_2d(gals.Colors, gals.Mag2s, gals.Mag2s, 'count', bins=2000)
    fig, ax = plt.subplots()
    #im = ax.imshow(np.log10(N.T), origin='lower',
    ##                           extent=[xedges[0], xedges[-1], yedges[0],
    #                                   yedges[-1]],
    #                           aspect='auto', interpolation='nearest',
    #                           cmap=plt.cm.RdBu)

    extent = [galshess[0][0], galshess[0][-1], galshess[1][-1], galshess[1][0]]
    imshow_kw={'norm': LogNorm(vmin=None, vmax=galshess[2].max()),
               'cmap': plt.cm.gray_r, 'interpolation': 'nearest',
               'aspect': 'equal', 'extent': extent}

    ax.plot(gals.Colorso, gals.Mag2o, ',', color='black')
    ax.autoscale(False)
    ax.set_xlim(-1, 3)
    ax.set_ylim(ylim)
    #ax.imshow(galshess[2], **imshow_kw)
    ax.set_aspect(1./2.)

    ax.set_xlabel('$%s-%s$' % (filter1, gals.galaxies[0].filter2), fontsize=20)
    ax.set_ylabel('$%s$' % gals.galaxies[0].filter2, fontsize=20)
    ax.hlines(mean_trgb, *ax.get_xlim(), lw=2, color='red', zorder=100)
    #ax.hlines(mean_trgb+1.5, *ax.get_xlim(), lw=2, color='red', zorder=100)
    #for k, v in poly_dict.items():
    #    ax.plot(v[:,0], v[:,1])
    #plt.colorbar(cs)
    #plt.savefig('opt_cmd_f606w.png', dpi=150)
    return fig, ax, gals


def plot_cum_sum_sfr(targets, file_origin='match-hmc'):
    '''cumulative sfr plot from match, no errors.'''
    match_sfh_src = snap_src + '/data/sfh_parsec/'

    fig, ax = plt.subplots()
    ngals = len(targets)
    bmap = brewer2mpl.get_map('Spectral', 'Diverging', ngals)
    cols = bmap.mpl_colors
    cols = ['#8ca8ba', '#0a6277', '#6a0c0c', '#bc741e', '#448833', '#88994b',
            '#89360f', '#b85121', '#aa4400']

    for i, target in enumerate(targets):
        match_sfh_file, = rsp.fileIO.get_files(match_sfh_src, '%s*sfh' % target.lower().replace('-deep', ''))
        sfh = sfh_tests.StarFormationHistories(match_sfh_file, file_origin=file_origin)
        age = 10**((sfh.data.lagef + sfh.data.lagei)/2. - 9)
        csfh = np.append(sfh.data.csfr, 0)
        ax.plot(age, csfh[1:], color=cols[i], lw=3,
                label='$%s$' % target.upper().replace('-','\!-\!'))

    #ax.set_xscale('log')
    ax.set_xlim(13.33, 0)
    ax.set_xlabel('$\\rm{Time\ (Gyr)}$', fontsize=20)
    ax.set_ylabel('$\\rm{Culmulative\ SF}$', fontsize=20)
    plt.legend(loc=0, frameon=False)
    plt.tick_params(labelsize=16)
    plt.savefig('csfr_ancients.png', dpi=150)


def plot_chi2_tests():
    gs =  gridspec.GridSpec(8, 1, width_ratios=[1, 3])
    for ax in gs:
        ax.set_xlabel('%s' % gal.target.upper().replace('-','\!-\!'))
        nmodels = len()
        ax.plot(1, )


def plot_cmd_lf(target, band):
    '''simple figure with the data and LF'''
    import model_plots
    lims = model_plots.load_plot_limits()
    row = lims[lims['target'] == target.lower()]
    target = target.upper()

    print row
    if band == 'opt':
        fits_src = snap_src + '/data/angst_no_trim'
        cmd_errors_kw = {}
        ymin, ymax = row['opt_cmdmin'], row['opt_cmdmax']
    else:
        fits_src = 'default'
        cmd_errors_kw = {'errclr': -.5}
        ymin, ymax = row['ir_cmdmin'], row['ir_cmdmax']
    fig = plt.figure(figsize=(6, 6))
    gs = gridspec.GridSpec(1, 2, width_ratios=[3, 1])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    gal = galaxy_tests.load_galaxy(target, band=band, fits_src=fits_src)

    gal.plot_cmd(gal.color, gal.mag2, ax=ax1, scatter_off=True)
    hist, bins = rsp.math_utils.hist_it_up(gal.mag2)
    ax1.set_ylim(ymax, ymin)
    if band == 'opt':
        xmin = -2
    if band == 'ir':
        xmin = -1

    #ax1.set_ylim()
    #ax2.set_ylim()

    ax1.set_xlim(xmin, ax1.get_xlim()[1])
    gal.decorate_cmd(ax=ax1, trgb=True, cmd_errors_kw=cmd_errors_kw)
    err = np.sqrt(hist)
    ax2.set_xscale('log')
    ax2.errorbar(hist, bins[:-1], xerr=err, color='black', drawstyle='steps-mid')
    ax2.plot(hist, bins[:-1], linestyle='mid-steps', color='black', lw=3)
    ax2.set_ylim(ax1.get_ylim())
    gal.put_a_line_on_it(ax2, gal.trgb)
    plt.subplots_adjust(wspace=0)
    ax2.tick_params(labelleft=False, labelright=True)
    ax2.set_xlim(min(hist), ax2.get_xlim()[1])
    ax2.set_xlabel('$\#$', fontsize=20)
    ax1.set_ylabel('$%s$' % gal.filter2, fontsize=20)
    ax1.set_xlabel('$%s-%s$' % (gal.filter1, gal.filter2), fontsize=20)
    add_color_cuts(gal.filter1, ax=ax1, vline_kw={'ls': '--', 'lw': 2})
    plt.tick_params(labelsize=16)
    outfname = '%s_%s_cmd.png' % (target, band)
    outfile = os.path.join(snap_src, 'plots', outfname)
    fig.savefig(outfile, dpi=150)

def add_color_cuts(filter1, ax=None, vline_kw=None):
    vline_kw = vline_kw or {}
    if ax is None:
        fig, ax = plt.subplots()
    color_cut = sfh_tests.get_color_cut(filter1)
    ax.vlines(color_cut, *ax.get_ylim(), **vline_kw)
    return ax

def add_completeness(target, ax=None, vline_kw=None):
    comp90 = sfh_tests.read_completeness_table()
    ind, = np.nonzero(comp90['target'] == target.upper())
    
    
if __name__ == '__main__':
    targets = galaxy_tests.load_targets('ancients')
    #plot_cum_sum_sfr(targets)
    [[plot_cmd_lf(target, band) for target in targets] for band in ['opt', 'ir']]


# below here was thesis spaz

def compare_mass_loss(mass=1.0, z=0.002, sets=['S_APR13', 'S_APR13VW93', 'S_MAR13'],
                      xcol='ageyr', ycol='dMdt'):
    '''
    made to plot a comparison between several mass prescriptions.
    Labels for the plot are set up stupidly, maybe in in_dict or labels arg...
    '''
    if xcol == 'ageyr':
        norm = 1e6
    else:
        norm = 1.                

    agb_tracks_dir = research_path + 'AGBTracks/CAF09'
    direcs = []
    labels = []
    for set in sets:
        if set == 'S_APR13':
            label = 'BS95'
        if set == 'S_APR13VW93':
            label = 'VW93'
        if set == 'S_MAR13':
            label = 'M13'
        direc = os.path.join(agb_tracks_dir, set)
        direc, = [os.path.join(direc,d) for d in os.listdir(direc) if str(z) in d]
        direc, = rsp.fileIO.get_files(direc, 'agb_%.2f*' % mass)
        direcs.append(direc)
        labels.append('$%s$' % label)
    tracks = [fileIO.get_numeric_data(t) for t in direcs]
    fig, ax = plt.subplots()
    for i in range(len(tracks)):
        ax.plot(tracks[i].data_array[xcol]/norm, tracks[i].data_array[ycol],
                label=labels[i], lw=2)
    ax.legend(loc=4, frameon=False)
    ax.set_xlabel('$Age\ (10^6 yr)$', fontsize=20)
    ax.set_ylabel('$\dot{M}$', fontsize=20)
    ax.text(.95, .90, '$M=%.2fM_\odot$' % mass, fontsize=16, transform=ax.transAxes, ha='right')
    ax.tick_params(labelsize=16)
    plt.savefig('compare_%s_loss_m%.2f.png' % (ycol, mass), dpi=150)
    return ax
    
    
def tp_initial_conditions():
    init_cond_base = research_path + 'AGBTracks/CAF09/S12_FIRST_TP/'
    search_term = '*.INP'
    init_conds = rsp.fileIO.get_files(init_cond_base, search_term)

    inps = [rsp.fileIO.readfile(ic, col_key_line=-1) for ic in init_conds]
    masses = [0.8, 1., 1.5, 2., 3., 4., 5.]
    z = 0.002
    
    fmt = ''
    line = ' & '.join(['%(m1).1f', '%(l1).3f', '%(te1).3f', '%(t1).2e'])
    line += ' \\\ \n'
    for inp in inps:
        if inp['z0'][0] != z:
            continue
        for mass in masses:
            ind, = np.nonzero(inp['m1'] == mass)
            fmt += (line % inp[ind])
    
    fmt.replace('e+0', '\times10^')
    header = r'''\begin{table}
\begin{tabular}{%s}
Mass/$\msun$) & $\log L/\Lsun$ & $\log T_{eff}$ (K) & Age (yr) \\
'''
    footer = r'''
\label{tab_init_cond}
\end{tabular}
\tablecomments{All data from \citet{Bressan2012} Tracks}
\end{table}
'''
    with open(research_path + 'communication/TPAGB Paper/init_cond_tab.tex', 'w') as f:
        f.write(header % ''.join('l' * 4))
        f.write(fmt)
        f.write(footer)

def not_sure_what_this_was_for():
    angst_tab = rsp.angst_tables.AngstTables()
    snap_sample = research_path + 'SNAP/data/galaxies'


    opt_snap = rsp.fileIO.get_files(snap_sample, '*trim.fits')
    ir_snap = rsp.fileIO.get_files(snap_sample, '*IR*.fits')


    for opt in opt_snap:
        if opt == research_path + 'SNAP/data/galaxies/10182_SN-NGC2403-PR_F606W_F814W.gst.trim.fits':
            continue

        gal = rsp.Galaxies.galaxy(opt, hla=False, filetype='fitstable')
        fig, ax = gal.plot_cmd(gal.color, gal.mag2, threshold=25, levels=3, scatter_args={'alpha':1})

        mtrgb_err = angst_tab.__getattribute__(gal.target.replace('-', '_'))['%s,%s' % (gal.filter1, gal.filter2)]['mTRGB_err']
        xlim = np.asarray(ax.get_xlim())
        xlim[0] -= 1
        xlim[1] += 1

        ax.fill_between(np.arange(*xlim), gal.trgb - mtrgb_err, 
                                         gal.trgb + mtrgb_err, 
                                         color='red', alpha=0.3)

        #annotate_fmt = '$m_{%s} = val \pm %.2f$' % (gal.filter2, mtrgb_err)
        #gal.put_a_line_on_it(ax, gal.trgb, ls='-', annotate_fmt=annotate_fmt.replace('val', '%.2f'))
        gal.put_a_line_on_it(ax, gal.trgb, ls='-', annotate=False)
        gal.text_on_cmd()
        #ax.set_ylim(ax.get_ylim()[0],22.5)

def ave_galaxy_cmd(targets=None, filter1=None, filter2=None):
    if targets is None:
        targets = galaxy_tests.targets_z002()

    if filter1 is None:
        filter1 = 'F160W'

    if filter2 is None:
        filter2 = 'F110W'

    fig, ax = galaxy_tests.multi_galaxy_hess(targets=targets, make_hess=True)

    

def bolometric_correction_plot(filter1=None, filter2=None, tp_mass=1., logg_val=0,
                               bc_file=None, ax=None):
    '''
    two trilegal simulations, one at Fe/H = 0, the other at Z = 0.002
    '''
    if filter1 is None:
        filter1 = 'F160W'
    if filter2 is None:
        filter2 = 'F814W'
    if bc_file is None:
        bc_file = '/Users/phil/research/padova_apps/photom/bc_odfnew/wfc3snap/bctab_p00.dat'
    
    # read in the bc file
    bc = rsp.fileIO.readfile(bc_file, col_key_line=1) 
    # don't include WD and M-giants
    bc = bc[:np.argmin(np.diff(bc['n'])) + 1]
    teff = bc['Teff']
    bc_mag1 = bc[filter1]
    bc_mag2 = bc[filter2]
    
    # read in the initial conditions file
    tp_inpfile = research_path + 'AGBTracks/CAF09/S12_FIRST_TP/S12_Z0.02_Y0.284_1TP.INP'    
    tp_inp = rsp.fileIO.readfile(tp_inpfile, col_key_line=-2) 

    # get the teff of the the inital mass tp_mass
    tp_masses = tp_inp['m1']
    tp_logtes = tp_inp['te1']
    tp_logls = tp_inp['l1']
    massind, = np.nonzero(tp_masses == tp_mass)
    tp_teff = 10 ** tp_logtes[massind]

    # select only relevant logg
    inds, = np.nonzero(bc['logg'] == logg_val)
    #loggs = np.unique(bc['logg'])
    #indss = [np.nonzero(bc['logg'] == logg) for logg in loggs]
    
    # make BC plot
    if ax is None:
        fig, ax = plt.subplots()
    #[ax.plot(teff[inds], -1 * bc_mag1[inds], label='$%s$' % filter1) for inds in indss]
    ax.plot(teff[inds], -1 * bc_mag1[inds], label='$%s$' % filter1, lw=2, color='k')
    ax.plot(teff[inds], -1 * bc_mag2[inds], '--', label='$%s$' % filter2, lw=2, color='k')
    #ax.vlines(tp_teff, *ax.get_ylim(), lw=2)
    ax.set_xlim(12000,3000)
    ax.set_ylim(1, -3)
    ax.legend(loc=0, frameon=False)
    ax.set_xlabel('$T_{eff}$', fontsize=20)
    ax.set_ylabel('$BC$', fontsize=20)
    #ax.set_xscale('log')
    return ax

def color_teff_plot(ax=None, filter_combos=None, logg_val=2, tracks_base=None,
                    prefix=None, photsys='wfc3snap', color_em=False):
    '''
    makes color vs teff plot with a line for each filter combo, if color_em
    is true, colors the tracks first.
    
    '''
    if color_em is True:
        pc.quick_color_em(track_base, prefix, photsys=photsys)

    if filter_combos is None:
        filter_combos = ['F475W-F814W',
                         'F555W-F814W',
                         'F606W-F814W',
                         'F110W-F160W']
    filters = np.unique(np.concatenate([f.split('-') for f in filter_combos]))

    search_term = '*F7_*PMS.dat.%s' % photsys    
    ts = pc.TrackSet(tracks_dir=tracks_base, prefix=prefix, track_search_term=search_term)
    
    # mass limit
    tinds = [i for i, m in enumerate(ts.masses) if (3. <= m <= 9.)]

    ts.squish(*np.concatenate([['LOG_TE'], filters]), **{'inds': tinds})

    # logg is not part of squish because not part of track data!
    loggs = np.concatenate([t.calc_logg() for t in np.asarray(ts.tracks)[tinds]])
    
    ginds, = np.nonzero(np.round(loggs, 0) == logg_val)
    
    teffs = 10 ** ts.LOG_TEs
    if ax is None:
        fig, ax = plt.subplots()
    
    for i, filts in enumerate(filter_combos):
        filt1, filt2 = filts.split('-')
        yval = ts.__getattribute__('%ss' % filt1) - ts.__getattribute__('%ss' % filt2)
        no_dupes, = np.nonzero(np.diff(yval) == 0)
        inds = list(set(no_dupes) & set(ginds))
        order = np.asarray(inds)[np.argsort(teffs[inds])]
        ax.plot(teffs[order], yval[order], lw=2, label='$%s$' % filts)

    ax.set_ylim(-.1, 4)
    ax.set_xlim(12000, 3000)
    ax.set_xlabel('$T_{eff}$', fontsize=20)
    ax.set_ylabel('$Color$', fontsize=20)
    ax.legend(loc=0, frameon=False)
    return ax

def two_panel():
    prefix = 'S12D_NS_Z0.02_Y0.284'
    tracks_base = '/Users/phil/research/parsec2match/S12_set/CAF09_S12D_NS'
    logg_val = 2
    fig, axs = plt.subplots(nrows=2)
    bolometric_correction_plot(tp_mass=1., logg_val=logg_val, ax=axs[0])
    color_teff_plot(tracks_base=tracks_base, prefix=prefix, ax=axs[1], logg_val=logg_val)
    plt.savefig('bc_plot.png', dpi=150)


def make_table(targets=None, deluxe=True):
    angst_data = rsp.angst_tables.AngstTables()
    tab2 = angst_data.snap_tab2
    tab3 = angst_data.snap_tab3

    if targets is None:
        targets = np.concatenate([galaxy_tests.targets_z002(), galaxy_tests.gi10_overlap()])

    fmt = ''
    line = '%(target)s & %(filters)s & %(Av).3f & %(mean_color).2f & '
    line += '%(mTRGB_F160W).3f$\pm$%(mTRGB_F160W_err).3f & %(50_completeness_F160W).2f & '
    line += '%(dmod).3f & %(logoh).3f$\pm$%(logoherr).3f & %(nrgb)s & %(nagb)i & '
    line += '%(ratio)s \\\ \n'

    for target in targets:
        atarget = target.replace('C-', 'C').replace('-', '_')
        datum = angst_data.__getattribute__(atarget)
        try:
            tab3_row, = tab3[np.nonzero(tab3['target'] == target)]        
            tab2_row, = tab2[np.nonzero(tab2['target'] == target)]
        except:
            btarget = target.replace('C-', 'C')
            tab3_row, = tab3[np.nonzero(tab3['target'] == btarget)]        
            tab2_row, = tab2[np.nonzero(tab2['target'] == btarget)]
        # turn the arrays into dicts
        tab3_dict = dict(zip(tab3.dtype.names, tab3_row.tolist()))
        tab2_dict = dict(zip(tab2.dtype.names, tab2_row.tolist()))
        # combine the dicts
        datum = dict((tab3_dict.items() + tab2_dict.items() + datum.items()))
        
        # extra stuff I want in the table
        filters = [key for key in datum.keys() if ',' in key]
        datum['filters'] = ','.join(np.unique(','.join(filters).split(',')))
        datum['Av'] = datum[filters[0]]['Av'] 
        datum['logoh'] = galaxy_tests.get_key_fromtable(target, 'logOH')
        datum['logoherr'] = galaxy_tests.get_key_fromtable(target, 'OHerr')
        #datum['Z'] = galaxy_tests.get_key_fromtable(target, 'Z')
        datum['nagb'] = galaxy_tests.get_key_fromtable(target, 'N_AGB')
        datum['nrgb'] = '...'
        datum['ratio'] = '...'
        fmt += line % datum
        
    fmt = fmt.replace('nan', '...')
    plain_header =  r'''\begin{table}
\begin{tabular}{%s}
Target/name & Opt.~Filters & $A_V$ & Mean NIR & $m_{\rm F160W}$ & $m_{\rm F160W}$ & $(m-M)_0$ & $\log [O/H]$ & $N_{\rm RGB}$ & $N_{\rm AGB}$ & $\frac{N_{\rm AGB}}{N_{\rm RGB}}$ \\
& & & Color & TRGB & \50%s complete & & & & & \\
1 & 2 & 3 & 4 & 5 & 6 & 7 & 8 & 9 & 10 & 11 \\
\hline
'''
    plain_footer = '''
\label{tab_sample}
\end{tabular}
\tablecomments{Columns 3--5 from \citet{Dalcanton2012}, column 8 from \citet[][and refs. therein]{Berg2012, Marble2010}}
\end{table}
'''

    deluxe_header = r'''\begin{deluxetable*}{%s}
\tablewidth{0pc}
\tabletypesize{\scriptsize}
\tablecaption{Galaxy Sample Parameters}
\tablehead{
    \colhead{Target/name} &
    \colhead{Opt.~Filters} &
    \colhead{$A_V$} &
    \colhead{Mean NIR} &
    \colhead{$m_{\rm F160W}$} &
    \colhead{$m_{\rm F160W}$ 50\%s} &
    \colhead{$(m-M)_0$} &
    \colhead{$\log [O/H]$} &
    \colhead{$N_{\rm RGB}$} &
    \colhead{$N_{\rm AGB}$} &
    \colhead{$\frac{N_{\rm AGB}}{N_{\rm RGB}}$} \\
    \colhead{} &
    \colhead{} &
    \colhead{} &
    \colhead{Color} &
    \colhead{TRGB} &
    \colhead{complete} &
    \colhead{} &
    \colhead{} &
    \colhead{} &
    \colhead{} &
    \colhead{} \\
    \colhead{1} &
    \colhead{2} &
    \colhead{3} &
    \colhead{4} &
    \colhead{5} &
    \colhead{6} &
    \colhead{7} &
    \colhead{8} &
    \colhead{9} &
    \colhead{10} &
    \colhead{11}
}
\startdata
'''

    deluxe_footer = r'''\enddata
\tablecomments{Columns 3--5 from \citet{Dalcanton2012}, column 8 from \citet[][and refs. therein]{Berg2012, Marble2010}}
\label{tab_sample}
\end{deluxetable*}
'''
    if deluxe is True:
        table_file = research_path + 'communication/TPAGB Paper/rgbagb2_deluxetable.tex'
        header = deluxe_header
        footer = deluxe_footer
    else:
        table_file = research_path + 'communication/TPAGB Paper/rgbagb2_table.tex'
        header = plain_header
        footer = plain_footer

    with open(research_path + 'communication/TPAGB Paper/rgbagb2_table.tex', 'w') as f:
        f.write(header % ('l' * 13, '%'))
        f.write(fmt)
        f.write(footer)

