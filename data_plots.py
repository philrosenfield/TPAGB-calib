import ResolvedStellarPops as rsp


'''
This is a work in progress, which sorts of intro figures should I have for the paper??

'''
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
