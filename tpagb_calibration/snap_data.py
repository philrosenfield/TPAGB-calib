import ResolvedStellarPops as rsp
import numpy as np
import os
import pyfits
from sfhs.vary_sfh import VarySFHs
import time
from IPython import parallel
import matplotlib.pylab as plt

def plot_lf_files():
    import plotting as pl
    from astropy.io import fits

    def translate_model_name(model, small=False):
        if 'oct' in model.lower():
            name = 'R75'
        if 'nov13eta' in model.lower():
            name = '\eta=0'
        if 'nov13_' in model.lower():
            name = 'mSC05'
        if 'feb' in model.lower():
            name = 'FEB14'
        if small is True:
            new_model = r'$\dot M_{\rm pre\!-\!dust}^{\rm %s}$' % name
        else:
            new_model = r'$\dot{M}_{\rm pre\!-\!dust}\!=\!%s$' % name
        return new_model

    opt_lfs = ['/Users/phil/Downloads/caf09_s_nov13_ddo71_opt_lf.dat',
               '/Users/phil/Downloads/caf09_s_nov13eta0_ddo71_opt_lf.dat',
               '/Users/phil/Downloads/caf09_s_oct13_ddo71_opt_lf.dat']

    ir_lfs = ['/Users/phil/Downloads/caf09_s_nov13_ddo71_ir_lf.dat',
              '/Users/phil/Downloads/caf09_s_nov13eta0_ddo71_ir_lf.dat',
              '/Users/phil/Downloads/caf09_s_oct13_ddo71_ir_lf.dat']
    colors = ['navy', 'orange', 'darkred']

    target = 'DDO71'
    opt_trgb = 23.71
    ir_trgb = 22.14
    opt_limit = 26.11
    ir_limit = 23.51

    base = '/Users/phil/research/TP-AGBcalib/SNAP/data/galaxies/'
    opt_galname = os.path.join(base, '9884_DDO71_F606W_F814W.gst.trim.fits')
    ir_galname = os.path.join(base, '11719_DDO71_IR_F110W_F160W.gst.fits')

    opt_gal = fits.open(opt_galname)
    ir_gal = fits.open(ir_galname)

    opt_mag2 = opt_gal[1].data['MAG2_ACS']
    ir_mag2 = ir_gal[1].data['MAG2_IR']

    _, (axs) = plt.subplots(ncols=2, figsize=(12, 6))
    plt.subplots_adjust(right=0.98, left=0.07, wspace=0.1)
    plot_data = True

    for i in range(len(opt_lfs)):
        plt_kw = {'color': colors[i], 'label': translate_model_name(opt_lfs[i])}
        axs = pl.plotting.compare_to_gal(opt_lfs[i], ir_lfs[i], opt_mag2,
                                         ir_mag2, 'F814W', 'F160W', opt_limit,
                                         ir_limit, opt_trgb=opt_trgb,
                                         ir_trgb=ir_trgb, axs=axs,
                                         plt_kw=plt_kw, target=target,
                                         trgb_excludes=[0.1, 0.2],
                                         plot_data=plot_data)
        plot_data = False

    axs[0].set_ylim(1, axs[0].get_ylim()[1])
    axs[0].set_xlim(21, 28.6)
    axs[1].set_xlim(19, 25.2)
    axs[1].set_ylim(1, axs[1].get_ylim()[1])
    plt.savefig('%s_lfs.eps' % target.lower())
    return axs


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
    nsfhs = 1
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

    clients = parallel.Client()
    clients.block = False
    here = os.getcwd()
    clients[:].execute('cd ~/research/TP-AGBcalib/code/TPAGB-calib/tpagb_calibration')
    clients[:].execute('from sfhs.vary_sfh import VarySFHs')
    clients[:].execute('cd %s' % here)

    res = [clients[i].apply(caller, vsfhs[i], vary_kws[i],)
           for i in range(len(sfh_files))]

    while False in [r.ready() for r in res]:
        time.sleep(10)
    print 'done.'
    os.system('ipcluster stop')


def rgb_agb_regions(offset, trgb_exclude, trgb, mag):
    # define RGB regions
    if offset >= 2.:
        low = offset
    else:
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

def normalize_to_flux(mag, flux, srgb, sagb):
    norm = flux / 10**(-.4 * mag)
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
    #with open(fname, 'r') as f:
    #    col_keys = f.readline().strip().split()

    # the columns I want
    #usecols = [col_keys.index(c) for c in cols]
    sgal.data = np.genfromtxt(fname, names=cols).view(np.recarray)
    return sgal


def select_rgbagb(StarPop, field1, field2, color_lim, offset, trgb_exclude, trgb):
    StarPop.color_cut, = np.nonzero(StarPop.data[field1] - StarPop.data[field2] > color_lim)
    StarPop.mag = StarPop.data[field2][StarPop.color_cut]
    StarPop.rgb, StarPop.agb = rgb_agb_regions(offset, trgb_exclude, trgb, StarPop.mag)
    return StarPop


def compare_model_data():
    agb_table = '/home/rosenfield/research/TP-AGBcalib/tables/melbourne2012_tab.dat'
    agb_tab = rsp.fileio.readfile(agb_table, string_column=0)

    color_lim = 0.1
    offsets = [23.51, 23.63, 23.85, 23.26, 22.67, 24.27]
    trgb_exclude = 0.2
    targets = ['DDO71', 'DDO78', 'HS117', 'KKH37', 'NGC2976-DEEP', 'SCL-DE1']
    model_folder =  '/home/rosenfield/research/TP-AGBcalib/SNAP/models/agb_flux'
    model_fmt = 'output_%s_%s_???.dat'
    data_folder =  '/home/rosenfield/research/TP-AGBcalib/SNAP/data/galaxies/'
    data_fmt = '11719_%s_IR_F110W_F160W.gst.fits'
    agb_mod = 'cmd_input_caf09_s_nov13'
    models = {}
    gals = []
    print 'target fmagb/fdagb fmagb/fdagb fdagb/fdtot, nagb_model/nagb_data agb/rgb_data agb_rgb_model'
    for i, target in enumerate(targets):
        dfname, = rsp.fileio.get_files(data_folder, data_fmt % target)
        mfnames = rsp.fileio.get_files(model_folder, model_fmt % (target, agb_mod))

        trgb  = rsp.angst_tables.angst_data.get_snap_trgb_av_dmod(target)[0]

        gal = rsp.StarPop()
        gal.target = target
        gal.data = pyfits.getdata(dfname)
        gal = select_rgbagb(gal, 'MAG1_IR', 'MAG2_IR', color_lim, offsets[i],
                            trgb_exclude, trgb)
        gal.agbflux = np.sum(10 ** (-.4 * gal.mag[gal.agb]))
        mag2 = gal.data['MAG2_IR']

        gal.totflux = np.sum(10 ** (-.4 * mag2))
        #gal.normflux = np.sum(10 ** (-.4 * gal.mag))
        #gal.totflux =  rsp.fileio.item_from_row(agb_tab, 'Galaxy', target, 'F160W_Flux')
        #agb_flux = rsp.fileio.item_from_row(agb_tab, 'Galaxy', target, 'fAGB_ftot') * gal.totflux
        #print gal.agbflux / agb_flux,  gal.agbflux, agb_flux, totflux, gal.totflux,  gal.totflux /  totflux
        gals.append(gal)

        sgals = []
        for mfname in mfnames:
            sgal = load_model(mfname)
            sgal.target = target
            sgal = select_rgbagb(sgal, 'F110W', 'F160W', color_lim, offsets[i],
                                 trgb_exclude, trgb)
            _,  sgal.norm_rgb, sgal.norm_agb = \
                normalize_simulation(sgal.mag, float(len(gal.rgb)), sgal.rgb,
                                     sgal.agb)
            gal.frgb = np.sum(10 ** (-.4 * gal.mag[gal.rgb]))
            sgal.frgb = np.sum(10 ** (-.4 * sgal.mag[sgal.rgb]))
            #print gal.frgb, sgal.frgb, gal.frgb/sgal.frgb
            sgal.norm  = gal.frgb / sgal.frgb
            sgal.totflux = np.sum(10 ** (-.4 * sgal.data.F160W))
            #sgal.norm = gal.totflux / sgal.totflux
            sgal.agbflux1 = np.sum(10 ** (-.4 * sgal.mag[sgal.agb])) * sgal.norm
            sgal.agbflux = np.sum(10 ** (-.4 * sgal.mag[sgal.norm_agb])) #* sgal.norm

            sgals.append(sgal)
            del sgal

        mean_flux1 = np.mean([s.agbflux1 for s in sgals])
        mean_flux = np.mean([s.agbflux for s in sgals])

        mean_agb = np.mean([float(len(s.norm_agb)) for s in sgals])
        mean_rgb = np.mean([float(len(s.norm_rgb)) for s in sgals])
        ndata_agb = float(len(gal.agb))
        print '%s %.3g %.3g %.3g %.3g %.4g %.4g' % (target, mean_flux1 / gal.agbflux,
                                               mean_flux / gal.agbflux,
                                               gal.agbflux / gal.totflux,
                                               mean_agb / ndata_agb,
                                               ndata_agb / float(len(gal.rgb)),
                                               mean_agb / mean_rgb)
        del gal
        models[target + '_mean'] = mean_flux
        models[target + '_nagb'] = mean_agb

        models[target] = sgals
    return gals, models

if __name__ == '__main__':
    call_vsfh()
    compare_model_data()
