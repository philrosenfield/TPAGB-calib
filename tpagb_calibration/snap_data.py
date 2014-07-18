import ResolvedStellarPops as rsp
import numpy as np
from matplotlib.nxutils import points_inside_poly
import os
import pyfits
from sfhs.vary_sfh import VarySFHs
import time
from IPython import parallel

def caller(vsfh, vary_kws={}):
    return vsfh.vary_the_SFH(**vary_kws)


def call_vsfh():
    base = '/home/rosenfield/research/TP-AGBcalib/SNAP/'

    file_origin = 'match-hmc'
    sfh_files = [base + 'data/sfh_parsec/ddo71.phil.final.sfh',
                 base + 'data/sfh_parsec/ddo78.phil.final.sfh',
                 base + 'data/sfh_parsec/hs117.hmc.sfh',
                 base + 'data/sfh_parsec/kkh37.phil.final.sfh',
                 base + 'data/sfh_parsec/ngc2976.phil.final.sfh',
                 base + 'data/sfh_parsec/scl-de1.phil.final.sfh']

    # gal input template
    gal_input_templates = [base + 'input/input_DDO71.dat',
                           base + 'input/input_DDO78.dat',
                           base + 'input/input_HS117.dat',
                           base + 'input/input_KKH37.dat',
                           base + 'input/input_NGC2976-DEEP.dat',
                           base + 'input/input_SCL-DE1.dat']

    # number of times to run trilegal
    nsfhs = 10
    # agb_model to use
    cmd_input_file = 'cmd_input_CAF09_S_NOV13.dat'
    # columns to save from trilegal
    cols = ['F110W', 'F160W']

    vsfhs = []
    vary_kws = []
    for i in range(len(sfh_files)):
        vsfhs.append(VarySFHs(sfh_file=sfh_files[i],
                             galaxy_input=gal_input_templates[i],
                             target=gal_input_templates[i].replace('.dat', '').split('_')[-1],
                             nsfhs=nsfhs,
                             cmd_input_file=cmd_input_file, file_origin=file_origin))
        vary_kws.append({'do_norm': False, 'write_culled': True, 'cols': cols})
    #os.system('ipcluster start -n=%i' % len(sfh_files))
    #import time
    #time.sleep(30)

    clients = parallel.Client()
    clients.block = False
    here = os.getcwd()
    clients[:].execute('cd ~/research/TP-AGBcalib/code/TPAGB-calib/tpagb_calibration')
    clients[:].execute('from sfhs.vary_sfh import VarySFHs')
    clients[:].execute('cd %s' % here)

    res = [clients[i].apply(caller, vsfhs[i], vary_kws[i],)
           for i in range(len(sfh_files))]

    while False in [r.ready() for r in res]:
        print 'waiting 30s ...'
        time.sleep(30)
        print 'checking ...'
    print 'done.'
    os.system('ipcluster stop')

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

def load_model(fname, cols=['F110W', 'F160W']):
    sgal = rsp.StarPop()

    # all the columns
    with open(fname, 'r') as f:
        col_keys = f.readline().strip().split()

    # the columns I want
    usecols = [col_keys.index(c) for c in cols]
    sgal.data = np.genfromtxt(fname, usecols=usecols,
                              names=cols).view(np.recarray)
    return sgal

def write_model(fname, data):
    header = '# %s \n' % ' '.join(data.dtype.names)
    with open(fname, 'w') as f:
        f.write(header)
        np.savetxt(f, data, fmt='%.4f')


def compare_model_data():


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
