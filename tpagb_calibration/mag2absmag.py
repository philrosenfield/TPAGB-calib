"""
convert mag to abs mag and save to a new file.

Uses a table to find the target, filters, photometric file. Could all be done
with just the name of the hla file though...

Will need BCDIR set for the mag2Mag call to work.
"""
import numpy as np
import os
import ResolvedStellarPops as rsp

from astropy.io import fits
from ResolvedStellarPops.tpagb_path_config import tpagb_path

gal_file = os.path.join(tpagb_path, 'SNAP/tables/snap_galaxies.dat')

gal_table = rsp.fileio.readfile(gal_file, string_column=[0,-4, -3, -2,-1],
                                string_length=216)

for i in range(len(gal_table)):
    phot = gal_table[i]['opt_phot']
    target = gal_table[i]['target']
    filt1 = gal_table[i]['opt1']
    filt2 = gal_table[i]['opt2']
    photsys = phot.split('_')[3].replace('-', '_')

    new_fname = phot.replace('.fits', '_absmag.dat')

    hdu = fits.getdata(os.path.join(tpagb_path, 'SNAP/data/angst_no_trim/',
                                    gal_table[i]['opt_phot']))
    mag1 = hdu['MAG1_{}'.format(photsys.split('_')[0].upper())]
    mag2 = hdu['MAG2_{}'.format(photsys.split('_')[0].upper())]
    Mag1 = rsp.astronomy_utils.mag2Mag(mag1, filt1, photsys,
                                       target=gal_table[i]['target'].upper(),
                                       filter1=filt1, filter2=filt2)

    Mag2 = rsp.astronomy_utils.mag2Mag(mag2, filt2, photsys,
                                       target=gal_table[i]['target'],
                                       filter1=filt1, filter2=filt2)
    header = '# mag_{0} mag_{1} Mag_{0} Mag_{1} \n'.format(filt1, filt2)
    rsp.fileio.savetxt(new_fname, np.column_stack((mag1, mag2, Mag1, Mag2)),
                       loud=True)
