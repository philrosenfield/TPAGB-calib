from astropy import coordinates as coord
import ResolvedStellarPops as rsp
import ResolvedStellarPops.convertz as convertz
'''
actually for the LMC, you can simplify things a lot:

1- the distance and reddening are known, take them from Stefano's paper.
    you don't need to relocate the TRGB.

2- you don't need to renormalize the simulations to the RGB numbers, because
    that was already done by Stefano. Stefano provides the SFR in units of
    Msun/yr. Integrate this SFR(t) from 0 to 15 Gyr, and you have the total
    mass of stars ever formed in that region of the LMC. This total mass is
    given to TRILEGAL in the line:

1.0e6 10.0 # object_mass, object_dist: total mass inside field, distance

Provided that you use the same IMF and binary fraction as Stefano, you
shouldn't have any problem.

Doing so, the only data you actually need for the LMC is the 2MASS data
(which is 99% complete for the AGB stars), the SFR files Stefano gave you
already, and the values you find in Stefano's paper (distance, Av, etc).

Let me know in case this is not clear.
'''

def parse_stefano_sfr():    
    filename = '/Users/phil/research/TP-AGBcalib/LMC_Calib/SFR_LMC88.dat'
    outfile = '/Users/phil/research/TP-AGBcalib/LMC_Calib/tri_SFR_LMC88.dat'
    data = rsp.fileIO.readfile(filename, col_key_line=1)

    cen_age = data['Center_age_bin']
    
    sfr = data['med_SFR']

    z = convertz.convertz(feh=data['med_FeH'])[1]
    zmin = convertz.convertz(feh=data['min_FeH'])[1]
    zmax = convertz.convertz(feh=data['max_FeH'])[1]
    # even dispersions
    zdisp = [np.mean([(zmax[i] - z[i]),
                      (z[i] - zmin[i])]) for i in range(len(z))]

    fmt = '%.3f %.3f %.4f %.4f\n'
    np.savetxt(outfile, np.array([cen_age, sfr, z, zdisp]).T, fmt=fmt)
    
    object_mass = np.sum(data['med_MASS']) * 1e6
    return outfile, object_mass

def make_trilegal_sim(cmd_input=None):
    if cmd_input is None:
    cmd_input = '/Users/phil/research/TP-AGBcalib/cmd_inputfiles/cmd_input_CAF09_S_MAR13.dat'
    cmd_input = '/Users/phil/research/TP-AGBcalib/cmd_inputfiles/cmd_input_CAF09_S_APR13.dat'

    object_sfr_file, object_mass = parse_stefano_sfr()
    galaxy_input = object_sfr_file.replace('.dat','_galinp.dat')
    output = '%s_%s' % (object_sfr_file.replace('.dat',''), cmd_input.split('cmd_input_')[1])

    gal_inp_pars = {'file_mag': 'tab_mag_odfnew/tab_mag_wfc3snap.dat',
                    'mag_limit_val': 20.6,
                    'mag_num': 12,
                    'binary_kind': 1,
                    'binary_frac': 0.3,
                    'object_sfr_mult_factorA': 1e9,
                    'object_dist': 47643.,
                    'object_av': 0.2, 
                    'object_avkind': 0,
                    'object_mass': object_mass,
                    'object_sfr_file': object_sfr_file}
    gal_inp = rsp.fileIO.input_parameters(default_dict=rsp.TrilegalUtils.galaxy_input_dict())
    gal_inp.add_params(gal_inp_pars)
    gal_inp.write_params(galaxy_input, rsp.TrilegalUtils.galaxy_input_fmt())
    rsp.TrilegalUtils.run_trilegal(cmd_input, galaxy_input, output)


vmcagbs = '/Users/phil/research/TP-AGBcalib/LMC_Calib/VMCAGBS_fmt.dat'

def read_vmc_table(filename):
    dtype = [('DEJ2000d', '<f8'),
             ('RAJ2000d', '<f8'),
             ('num', '<f8'),
             ('seq', '<f8'),
             ('RAJ2000', '|S12'),
             ('DEJ2000', '|S12'),
             ('tau', '<f8'),
             ('logMdot', '<f8'),
             ('Lum', '<f8'),
             ('Cl', '|S4'),
             ('Q', '|S1'),
             ('P', '|S1'),
             ('SED', '|S3'),
             ('Umag', '<f8'),
             ('Bmag', '<f8'),
             ('Vmag', '<f8'),
             ('Imag', '<f8'),
             ('YmagV', '<f8'),
             ('JmagV', '<f8'),
             ('KsmagV', '<f8'),
             ('Jmag2', '<f8'),
             ('Hmag2', '<f8'),
             ('Ksmag2', '<f8'),
             ('u3.6_1', '<f8'),
             ('u3.6_2', '<f8'),
             ('u4.5_1', '<f8'),
             ('u4.5_2', '<f8'),
             ('u5.8_1', '<f8'),
             ('u5.8_2', '<f8'),
             ('u8.0_1', '<f8'),
             ('u8.0_2', '<f8'),
             ('u24_1', '<f8'),
             ('u24_2', '<f8'),
             ('chi2C', '<f8'),
             ('chi2O', '<f8')]
    return np.genfromtxt(filename, dtype=dtype)
    
def read_lmc_cat(filename):
    with open(filename, 'r') as f:
        header = f.readline()
    col_keys = header.replace('#', '').strip().split()
    
    return np.genfromtxt(filename, names=col_keys)

    
lmc_cat_name = '/Users/phil/research/TP-AGBcalib/LMC_Calib/LMC88_psf_Y_J_Ks.cat.update'
lmc_cat = read_lmc_cat(lmc_cat_name)

def make_plot(output, lmc_cat):
    sgal = rsp.Galaxies.simgalaxy(output, filter1='J', filter2='Ks')
    sgal.all_stages('TPAGB')
    fig, (ax1, ax2) = plt.subplots(ncols=2)
    ax1.plot(lmc_cat['J_MAG'] - lmc_cat['Ks_MAG'], lmc_cat['Ks_MAG'],',', color='black')
    ax2.plot(sgal.color, sgal.mag2, ',', color='black')
    ax2.plot(sgal.color[sgal.itpagb], sgal.mag2[sgal.itpagb], '.', color='red')
    ax2.set_ylim(ax2.get_ylim()[::-1])
    ax1.set_xlim(ax2.get_xlim())
    ax1.set_ylim(ax2.get_ylim())
    ax1.set_xlabel('$J-Ks$')
    ax1.set_ylabel('$Ks$')
    ax2.set_title('$N_{TPAGB} = %i$' % len(sgal.itpagb))
    [ax.tick_params(labelsize=16) for ax in [ax1, ax2]]
    plt.savefig(output.replace('.dat','.png'))


data = read_vmc_table(vmcagbs)
radec = [coord.FK5Coordinates(r, d) for r, d in zip(data['RAJ2000'], data['DEJ2000'])]
ra = np.array([r.ra.degrees for r in radec])
dec = np.array([r.dec.degrees for r in radec])

