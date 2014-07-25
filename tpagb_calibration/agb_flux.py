import ResolvedStellarPops as rsp
from matplotlib.nxutils import points_inside_poly
import os
import pyfits
%pylab

def rgb_agb_regions(offset, trgb_exclude, trgb, mag):
    # define RGB regions
    low = offset + trgb
    mid = trgb + trgb_exclude

    # Recovered stars in simulated RGB region.
    srgb = rsp.starpop.stars_in_region(mag, low, mid)

    # define AGB regions
    mmid = trgb - trgb_exclude
    high = 10.

    # Recovered stars in simulated AGB region.
    sagb = rsp.starpop.stars_in_region(mag, mmid, high)
    return srgb, sagb

def normalize_simulation(mag, ndata_rgb, srgb, sagb):
    norm = ndata_rgb / float(len(srgb))
    #print 'Normalization: %f' % norm

    # random sample the data distribution
    rands = np.random.random(len(mag))
    ind, = np.nonzero(rands < norm)

    # scaled rgb: norm + in rgb
    rgb = list(set(ind) & set(srgb))

    # scaled agb
    agb = list(set(ind) & set(sagb))
    return norm, rgb, agb

def load_model(fname):
    sgal = rsp.StarPop()

    # all the columns
    with open(fname, 'r') as f:
        col_keys = f.readline().strip().split()

    # the columns I want
    cols = ['F110W', 'F160W']
    usecols = [col_keys.index(c) for c in cols]
    sgal.data = np.genfromtxt(fname, usecols=usecols,
                              names=cols).view(np.recarray)
    return sgal

def select_rgbagb(StarPop, field1, field2, color_lim, offset, trgb_exclude, trgb):
    StarPop.color_cut, = np.nonzero(StarPop.data[field1] - StarPop.data[field2] > color_lim)
    StarPop.mag = StarPop.data[field2][StarPop.color_cut]
    StarPop.rgb, StarPop.agb = rgb_agb_regions(offset, trgb_exclude, trgb, StarPop.mag)
    return StarPop

offset = 1.5
trgb_exclude = 0.2
color_lim = 0.3

data_folder = '/home/rosenfield/research/TP-AGBcalib/SNAP/data/galaxies/'

data_files = ['11719_DDO71_IR_F110W_F160W.gst.fits',
              '11719_DDO78_IR_F110W_F160W.gst.fits',
              '11719_HS117_IR_F110W_F160W.gst.fits',
              '11719_KKH37_IR_F110W_F160W.gst.fits',
              '11719_NGC2976-DEEP_IR_F110W_F160W.gst.fits',
              '11719_SCL-DE1_IR_F110W_F160W.gst.fits']

model_folder =  '/home/rosenfield/research/TP-AGBcalib/SNAP/models/varysfh/comp_corr/'
model_files = ['ddo71/caf09_s_nov13/mc/output_DDO71_ir_lf.dat',
               'ddo78/caf09_s_nov13/mc/output_DDO78_029_opt_best.dat',
               'hs117/caf09_s_nov13/mc/output_HS117_027_opt_best.dat',
               'kkh37/caf09_s_nov13/mc/output_KKH37_027_opt_best.dat',
               'ngc2976-deep/caf09_s_nov13/mc/output_NGC2976-DEEP_036_opt_best.dat',
               'scl-de1/caf09_s_nov13/mc/output_SCL-DE1_015_opt_best.dat']

sgals = []
gals = []
for i in range(len(data_files)):
    dfname = os.path.join(data_folder, data_files[i])
    mfname =  os.path.join(model_folder, model_files[i])

    sgal = load_model(mfname)
    sgal.target = os.path.split(dfname)[1].split('_')[1]
    trgb  = rsp.angst_tables.angst_data.get_snap_trgb_av_dmod(sgal.target)[0]
    sgal = select_rgbagb(sgal, 'F110W', 'F160W', color_lim, offset, trgb_exclude, trgb)

    gal = rsp.StarPop()
    gal.target = sgal.target
    gal.data = pyfits.getdata(dfname)
    gal = select_rgbagb(gal, 'MAG1_IR', 'MAG2_IR', color_lim, offset, trgb_exclude, trgb)
    gal.agbflux = np.sum(10 ** (-.4 * gal.mag[gal.agb]))

    sgal.norm,  sgal.norm_rgb, sgal.norm_agb = normalize_simulation(sgal.mag, float(len(gal.rgb)), sgal.rgb, sgal.agb)
    sgal.agbflux = np.sum(10 ** (-.4 * sgal.mag[sgal.norm_agb]))
    print float(len(gal.rgb)) / float(len(sgal.norm_rgb))
    print '%s %.3g %i %i %.3g' % (sgal.target, sgal.agbflux/gal.agbflux, len(gal.agb), len(sgal.agb), float(len(sgal.agb))/float(len(gal.agb)))
    sgals.append(sgal)
    gals.append(gal)
