'''
This code will add a column "field" to each fits file (saving new ones)
with -1 to mean the data are outside the overlap regions, 0 meaning they
are in the overlap region, 6 if in field 6, 12 if in field 12, 15
if in field 15, and 18 if in both field 6 and field 12.
'''
import matplotlib.nxutils as nx
import matplotlib.pylab as plt
import numpy as np
import os
import pyfits

# global locations of data files.
home = os.environ['HOME']
data_src = os.path.join(home, 'research/TP-AGBcalib/PHAT/data/')

def get_verts(x, y, dx=None, dy=None, nbinsx=100, nbinsy=100, smooth=False):
    ''' get the approx. footprint of some 2d array '''
    ymin = y.min()
    ymax = y.max()
    xmin = x.min()
    xmax = x.max()

    if dx is None or dy is None:
        dx = (xmax - xmin) / nbinsx
        dy = (ymax - ymin) / nbinsy
    else:
        nbinsx = (xmax - xmin) / dx
        nbinsy = (ymax - ymin) / dy

    ymid = []
    min_x = []
    max_x = []
    for j in range(nbinsy):
        yinner = ymin + j * dy
        youter = ymin + (j + 1) * dy
        # counterintuitive because I'm dealing with mags...
        ind, = np.nonzero((y > yinner) & (y < youter))
        if len(ind) > 0:
            if smooth:
                min_x.append(np.average(x[ind]) - 3. * np.std(x[ind]))
                max_x.append(np.average(x[ind]) + 3. * np.std(x[ind]))
                ymid.append((yinner + youter) / 2.)
            else:
                min_x.append(np.min(x[ind]))
                max_x.append(np.max(x[ind]))
                ymid.append((yinner + youter) / 2.)

    max_x.reverse()
    ymidr = ymid[:]
    ymidr.reverse()

    # close polygon
    max_x.append(min_x[0])

    # close polygon
    ymidr.append(ymid[0])

    # get verticies of polygon
    xs = np.concatenate((min_x, max_x))
    ys = np.concatenate((ymid, ymidr))
    verts = np.column_stack((xs, ys))

    return verts


def add_col(fits, new_arr, column_name='FIELD', col_format='D'):
    '''pyfits can't add a frickin column. This makes a new table with the new_arr as a new column'''
    cols = [c for c in fits[1].columns]
    col = pyfits.Column(name=column_name, format=col_format, array=new_arr)
    cols.append(col)
    cols = pyfits.ColDefs(cols)
    return pyfits.new_table(cols)


def test_outputs(new_phat_data_sfhs, new_matched_name, martha=True):
    # testing
    for i in range(len(new_phat_data_sfhs)):
        phat_hdu = pyfits.open(new_phat_data_sfhs[i])
        phat_ra = phat_hdu[1].data.field('RA')
        phat_dec = phat_hdu[1].data.field('DEC')
        phat_field = phat_hdu[1].data.field('field')
        inds, = np.nonzero(phat_field>0)
        plt.plot(phat_ra[inds], phat_dec[inds], ',')

    matched_hdu = pyfits.open(new_matched_name)
    if martha:
        ra = 'PHAT_RA'
        dec = 'PHAT_DEC'
    else:
        ra = 'RA'
        dec = 'DEC'
    matched_ra = matched_hdu[1].data.field(ra)
    matched_dec = matched_hdu[1].data.field(dec)
    matched_field = matched_hdu[1].data.field('FIELD')
    fields = [6., 12., 15., 18.]  # field=18 is because field 6 and 12 overlap.
    indss = [np.nonzero(matched_field == 6.),
             np.nonzero(matched_field == 12.),
             np.nonzero(matched_field == 15.),
             np.nonzero(matched_field == 18.)]
    for i in range(len(indss)):
        plt.plot(matched_ra[indss[i]], matched_dec[indss[i]], '.', label=fields[i])
    plt.legend()
    plt.show()


def main(matched_name='defualt', martha=True):
    # orignial phat fields that match ran SFHs...
    phat_data_sfhs = [os.path.join(data_src, '12055_M31-B21-F06.gst.fits'),
                      os.path.join(data_src, '12055_M31-B21-F12.gst.fits'),
                      os.path.join(data_src, '12055_M31-B21-F15.gst.fits')]

    # phat and irac fits file from Martha
    martha = True # only to use PHAT_RA PHAT_DEC instead of RA DEC field names
    matched_name = home + '/research/TP-AGBcalib/PHAT/data/Spitzer/B21_phat+irac_v2.fits'

    #martha = False
    #matched_name = home + '/research/TP-AGBcalib/PHAT/data/Spitzer/b21-6filt-cut-shallow.fits'

    # name of the new phat-irac data (with added "field" column)
    new_matched_name = matched_name.replace('.fits', '_fields.fits').replace('+', '_')

    # names of the new phat data (with added "field" column)
    new_phat_data_sfhs = [d.replace('.gst.fits', '_trimmed_ir.gst.fits') for d in phat_data_sfhs]

    # fields will be floats labeled to each star in the fits files
    fields = [int(os.path.split(f)[1].split('-')[-1].split('.')[0].replace('F', '')) for f in phat_data_sfhs]

    # load phat-irac data
    matched_hdu = pyfits.open(matched_name)
    if martha:
      ra = 'PHAT_RA'
      dec = 'PHAT_DEC'
    else:
      ra = 'RA'
      dec = 'DEC'

    matched_ra = matched_hdu[1].data.field(ra)
    matched_dec = matched_hdu[1].data.field(dec)

    # unmatched gets some huge - number RA/Dec
    inds, = np.nonzero((matched_ra > 0) & (matched_dec > 0))

    # new array will assign stars out of the regions as field=-1.
    new_marr = np.zeros(len(matched_dec)) - 1.
    # matched inds will have field = 0
    new_marr[inds] += 1.

    matched_ra = matched_ra[inds]
    matched_dec = matched_dec[inds]

    matched_verts = get_verts(matched_ra, matched_dec, nbinsx=100, nbinsy=20)
    matched_data = np.column_stack([matched_ra, matched_dec])

    for i in range(len(phat_data_sfhs)):
      # load phat data
      phat_hdu = pyfits.open(phat_data_sfhs[i])
      phat_ra = phat_hdu[1].data.field('RA')
      phat_dec = phat_hdu[1].data.field('DEC')

      # new array will assign stars out of the regions as field=-1.
      new_parr = np.zeros(len(phat_dec)) - 1.

      phat_verts = get_verts(phat_ra, phat_dec, nbinsx=100, nbinsy=100)
      phat_data = np.column_stack([phat_ra, phat_dec])

      # phat data that overlaps with matched_verts
      pinds, = np.nonzero(nx.points_inside_poly(phat_data, matched_verts))
      new_parr[pinds] += fields[i] + 1  # new_parr is initialized at - 1
      new_phat_fits = add_col(phat_hdu, new_parr)

      new_phat_fits.writeto(new_phat_data_sfhs[i], clobber=True)
      # matched data that overlaps with phat_verts
      minds, = np.nonzero(nx.points_inside_poly(matched_data, phat_verts))
      new_marr[inds[minds]] += fields[i]

    new_matched_fits = add_col(matched_hdu, new_marr)
    new_matched_fits.writeto(new_matched_name, clobber=True)

    # plot the new files for visual inspection
    test_outputs(new_phat_data_sfhs, new_matched_name, martha=martha)

if __name__ == '__main__':
    # could add argvs here ...
    main()