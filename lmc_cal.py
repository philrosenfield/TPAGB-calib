from astropy import coordinates as coord

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

data = read_vmc_table(vmcagbs)
radec = [coord.FK5Coordinates(r, d) for r, d in zip(data['RAJ2000'], data['DEJ2000'])]
ra = np.array([r.ra.degrees for r in radec])
dec = np.array([r.dec.degrees for r in radec])


