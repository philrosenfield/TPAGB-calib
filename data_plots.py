import ResolvedStellarPops as rsp


'''
This is a work in progress, which sorts of intro figures should I have for the paper??

'''
def not_sure_what_this_was_for():
    angst_tab = rsp.angst_tables.AngstTables()
    snap_sample = '/Users/phil/research/TP-AGBcalib/SNAP/data/galaxies'


    opt_snap = rsp.fileIO.get_files(snap_sample, '*trim.fits')
    ir_snap = rsp.fileIO.get_files(snap_sample, '*IR*.fits')


    for opt in opt_snap:
        if opt == '/Users/phil/research/TP-AGBcalib/SNAP/data/galaxies/10182_SN-NGC2403-PR_F606W_F814W.gst.trim.fits':
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

def bolometric_correction_plot(sim_file=None, filter1=None, filter2=None,
                               tp_mass=1., logg=0, bc_file=None):
    '''
    two trilegal simulations, one at Fe/H = 0, the other at Z = 0.002
    '''
    if filter1 is None:
        filter1 = 'F160W'
    if filter2 is None:
        filter2 = 'F814W'
    if bc_file is None:
        bc_file = '/Users/phil/research/padova_apps/photom/bc_odfnew/wfc3snap/bctab_p00.dat'
    
    bc = rsp.fileIO.readfile(bc_file, col_key_line=1) 
    # don't include WD and M-giants
    bc = bc[:np.argmin(np.diff(bc['n'])) + 1]
    tp_inpfile = '/Users/phil/research/TP-AGBcalib/AGBTracks/CAF09/S12_FIRST_TP/S12_Z0.02_Y0.284_1TP.INP'    
    tp_inp = rsp.fileIO.readfile(tp_inpfile, col_key_line=-2) 

    tp_masses = tp_inp['m1']
    tp_logtes = tp_inp['te1']
    tp_logls = tp_inp['l1']
    massind, = np.nonzero(tp_masses == tp_mass)
    tp_teff = 10 ** tp_logtes[massind]

    #inds, = np.nonzero(bc['logg'] == logg)
    loggs = np.unique(bc['logg'])
    indss = [np.nonzero(bc['logg'] == logg) for logg in loggs]
    
    teff = bc['Teff']
    bc_mag1 = bc[filter1]
    #bc_mag2 = bc[filter2]

    fig, ax = plt.subplots()
    [ax.plot(teff[inds], -1 * bc_mag1[inds], label='$%s$' % filter1) for inds in indss]
    #ax.plot(teff[inds], -1 * bc_mag1[inds], label='$%s$' % filter1, lw=2, color='k')
    #ax.plot(teff[inds], -1 * bc_mag2[inds], '--', label='$%s$' % filter2, lw=2, color='k')
    #ax.vlines(tp_teff, *ax.get_ylim(), lw=2)
    ax.set_xlim(12000,3000)
    ax.set_ylim(1,-3)
    #ax.set_xscale('log')
