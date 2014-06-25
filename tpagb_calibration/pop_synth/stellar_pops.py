import numpy as np

def rgb_agb_regions(sgal, offsets, trgb_excludes, opt_trgb,
                    ir_trgb, opt_mag, ir_mag):
    # define RGB regions
    opt_low = offsets[0]
    opt_mid = opt_trgb + trgb_excludes[0]

    ir_low = offsets[1]
    ir_mid = ir_trgb + trgb_excludes[1]

    # Recovered stars in simulated RGB region.
    sopt_rgb = sgal.stars_in_region(opt_mag, opt_low, opt_mid)
    sir_rgb = sgal.stars_in_region(ir_mag, ir_low, ir_mid)

    # define AGB regions
    opt_mid = opt_trgb - trgb_excludes[0]
    opt_high = 10.

    ir_mid = ir_trgb - trgb_excludes[1]
    ir_high = 10.

    # Recovered stars in simulated AGB region.
    sopt_agb = sgal.stars_in_region(opt_mag, opt_mid, opt_high)
    sir_agb = sgal.stars_in_region(ir_mag, ir_mid, ir_high)
    return sopt_rgb, sir_rgb, sopt_agb, sir_agb


def normalize_simulation(opt_mag, ir_mag, nopt_rgb, nir_rgb, sopt_rgb, sir_rgb,
                         sopt_agb, sir_agb):
    opt_norm = nopt_rgb / float(len(sopt_rgb))
    ir_norm = nir_rgb / float(len(sir_rgb))

    print 'OPT Normalization: %f' % opt_norm
    print 'IR Normalization: %f' % ir_norm

    # random sample the data distribution
    rands = np.random.random(len(opt_mag))
    opt_ind, = np.nonzero(rands < opt_norm)
    rands = np.random.random(len(ir_mag))
    ir_ind, = np.nonzero(rands < ir_norm)

    # scaled rgb: norm + in rgb
    opt_rgb = list(set(opt_ind) & set(sopt_rgb))
    ir_rgb = list(set(ir_ind) & set(sir_rgb))

    # scaled agb
    opt_agb = list(set(opt_ind) & set(sopt_agb))
    ir_agb = list(set(ir_ind) & set(sir_agb))
    return opt_norm, ir_norm, opt_rgb, ir_rgb, opt_agb, ir_agb
