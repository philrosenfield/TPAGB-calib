import itertools
from TPAGBparams import phat_src
import os
import ResolvedStellarPops as rsp
import numpy as np
from sfhs.vary_sfh import VarySFHs
from pop_synth.stellar_pops import rgb_agb_regions
import calibrate
from plotting import Plotting
import time



__all__ = ['prepare_data', 'prepare_vsfh_run']

def load_data_files():
    base = '/home/rosenfield/research/TP-AGBcalib/PHAT/data/matchphot'
    opt_names = ['M31-B21_3x6-006.gst.match',
                 'M31-B21_3x6-012.gst.match',
                 'M31-B21_3x6-015.gst.match']
    opt_files = [os.path.join(base, opt_name) for opt_name in opt_names]

    ir_names = ['12055_M31-B21-F06-IR_F110W_F160W.gst.match',
                '12055_M31-B21-F12-IR_F110W_F160W.gst.match',
                '12055_M31-B21-F15-IR_F110W_F160W.gst.match']
    ir_files = [os.path.join(base, ir_name) for ir_name in ir_names]
    return opt_files, ir_files

def prepare_data_table(offsets=[-2., -1.5], trgb_excludes=[.1, .2]):
    opt_color_min = 1.5
    ir_color_min = 0.6
    Mtrgb = -4.05
    dmod = 24.47
    Av = 0.4
    photsys = 'phat_agb'
    ir_trgb = rsp.astronomy_utils.Mag2mag(Mtrgb, 'F160W', photsys,
                                                  dmod=dmod, Av=Av)
    opt_trgb = rsp.astronomy_utils.Mag2mag(Mtrgb, 'F814W', photsys,
                                                  dmod=dmod, Av=Av)
    offsets[0] = opt_trgb - offsets[0]
    offsets[1] = ir_trgb - offsets[1]

    fmt = '%s %i %i %i %i %.2f %.2f %.2f %.2f \n'
    header = '# opt: F814W ir: F160W Mtrgb: %.2f dmod: %.2f Av: %.1f \n'
    header += '# opt, ir: color_min (%.1f, %.1f) trgb (%.2f, %.2f) '
    header += 'offsets (%.1f, %.1f) trgb exclude (%.1f, %.1f) \n'
    header += '# target nopt_rgb nopt_agb nir_rgb nir_agb opt_max ir_max'
    header += ' opt_min ir_min \n'

    out = open('phat_b21_data.dat', 'w')
    out.write(header % (Mtrgb, dmod, Av, opt_color_min, ir_color_min,  opt_trgb,
                        ir_trgb, offsets[0], offsets[1], trgb_excludes[0],
                        trgb_excludes[1]))

    opt_files, ir_files = load_data_files()
    for i in range(len(opt_files)):
        opt_gal = rsp.galaxies.galaxy.Galaxy(opt_files[i],
                                         filetype='match_phot', angst=False,
                                         hla=False, filter1='F475W',
                                         filter2='F814W')
        ir_gal = rsp.galaxies.galaxy.Galaxy(ir_files[i],
                                         filetype='match_phot', angst=False,
                                         hla=False, filter1='F110W',
                                         filter2='F160W')
        opt_color_cut, = np.nonzero((opt_gal.color) > opt_color_min)
        ir_color_cut, = np.nonzero((ir_gal.color) > ir_color_min)
        ir_mag = ir_gal.mag2[ir_color_cut]
        opt_mag = opt_gal.mag2[opt_color_cut]
        opt_rgb, ir_rgb, opt_agb, ir_agb = rgb_agb_regions(opt_gal, offsets,
                                                          trgb_excludes,
                                                          opt_trgb,
                                                          ir_trgb,
                                                          opt_mag,
                                                          ir_mag)
        out.write(fmt % (opt_files[i].split('.')[0], len(opt_rgb), len(opt_agb),
                         len(ir_rgb), len(ir_agb), opt_gal.mag2.max(),
                         ir_gal.mag2.max(), opt_gal.mag2.min(),
                         ir_gal.mag2.min()))


def prepare_vsfh_run(nsfhs):
    '''
    Run a number of SFH variations on a galaxy.
    If passed default to the args, will attempt to find the file based on the
    galaxy name.

    ARGS:
    galaxy_name: target name, ex: ddo71 (case doesn't matter)
    cmd_input_file: filename, ex: 'cmd_input_CAF09_S_OCT13.dat'
    match_sfh_file: 'default', the sfh file from match.
    match_fileorigin: 'match-grid', which type of sfh file
                        'match' or 'match-grid'
    galaxy_input_file: 'default', base input file for trilegal, will be copied.

    mk_tri_sfh_kw: A dict to be passed to VarySFH.make_trilegal_sfh
        default: random_sfh = True, random_z = False

    make_many_kw: A dict to be passed to VarySFH.prepare_trilegal_sfr
        default: nsfhs = 50, mk_tri_sfh_kw dict.

    vary_sfh_kw: A dict to be passed to VarySFH.vary_the_sfh
        default: diag_plots = True, make_many_kw dict

    RETURNS:
    VarySFHs class
    '''
    track_model = 'PARSEC'
    base = phat_src

    # Data information
    targets = ['M31-B21_3x6-006',
               'M31-B21_3x6-012',
               'M31-B21_3x6-015']

    dmod = 24.47
    Av = 0.
    photsys = 'phat_agb'
    opt_color_min = 1.5
    ir_color_min = 0.6

    # Match SFH files
    sfh_src = os.path.join(base, 'match', 'B21_AGB_%s' % track_model)
    sfh_ext = '.sfh.mcmc.zc'

    # galaxy input file for trilegal
    galaxy_input_scr = os.path.join(base, 'vary_sfh', 'input')
    galaxy_input_ext = '.inp'  # could also be .aringer.inp

    cmd_input_files = ['cmd_input_CAF09_S_NOV13.dat']

    table_file = os.path.join(base, 'tables', 'phat_b21_data.dat')
    dtype=[('target', '|S25'), ('nopt_rgb', '<f8'), ('nopt_agb', '<f8'),
           ('nir_rgb', '<f8'), ('nir_agb', '<f8'), ('opt_max', '<f8'),
           ('ir_max', '<f8'), ('opt_min', '<f8'), ('ir_min', '<f8')]
    m31_data = np.genfromtxt(table_file, dtype=dtype)
    dmag = 0.1

    vsfh_kw = {'file_origin': 'match-hmc',
               'nsfhs': nsfhs,
               'filter1': 'F475W',
               'photsys': photsys,
               'opt_color_min': opt_color_min,
               'ir_color_min': ir_color_min,
               'dmod': dmod,
               'Av': Av,
               'Mtrgb': -4.05}
    vsfh_kws = []
    vsfhs = []

    for target, cmd_input in itertools.product(targets, cmd_input_files):
        agb_mod = cmd_input.replace('cmd_input_', '').lower().split('.')[0]

        galaxy_input = os.path.join(galaxy_input_scr, target + galaxy_input_ext)
        sfh_file = os.path.join(sfh_src, target + sfh_ext)

        outfile_loc = os.path.join(base, 'vary_sfh', target, agb_mod)
        rsp.fileIO.ensure_dir(outfile_loc)

        data_line = m31_data['target'==target]
        opt_bins = np.arange(data_line['opt_min'], data_line['opt_max'], dmag)
        ir_bins = np.arange(data_line['ir_min'], data_line['ir_max'], dmag)

        #print target, cmd_input
        vsfh_kw.update({'target': target.lower(),
                        'galaxy_input': galaxy_input,
                        'sfh_file': sfh_file,
                        'outfile_loc': outfile_loc,
                        'cmd_input_file': cmd_input,
                        'nopt_rgb': data_line['nopt_rgb'],
                        'nir_rgb': data_line['nir_rgb'],
                        'nopt_agb': data_line['nopt_agb'],
                        'nir_agb': data_line['nir_agb'],
                        'opt_bins': opt_bins,
                        'ir_bins': ir_bins})

        vsfhs.append(VarySFHs(**vsfh_kw))
        vsfh_kws.append(vsfh_kw)
    return vsfhs, vsfh_kws

def phat_data(nsfhs=50, nprocs=6):
    vsfhs, vsfh_kws = prepare_vsfh_run(nsfhs)
    nruns = len(vsfhs) * nsfhs
    if nruns > 10:
        os.system('ipcluster start -n=%i' % nprocs)
        time.wait(45)
        calibrate.main(vsfhs, vsfh_kws=vsfh_kws)
    os.system('ipcluster stop')
    #results = calibrate.main(vsfhs, vsfh_kws=vsfh_kws)
    opt_files, ir_files = load_data_files()
    for i, vsfh in enumerate(vsfhs):
        if nruns < 10:
            vsfh.vary_the_SFH()
            vsfh.write_results()
        pl = Plotting(vsfh)
        opt_gal = rsp.galaxies.galaxy.Galaxy(opt_files[i],
                                         filetype='match_phot', angst=False,
                                         hla=False, filter1='F475W',
                                         filter2='F814W')
        ir_gal = rsp.galaxies.galaxy.Galaxy(ir_files[i],
                                         filetype='match_phot', angst=False,
                                         hla=False, filter1='F110W',
                                         filter2='F160W')

        pl.compare_to_gal(opt_gal, ir_gal, 28, 28, narratio=False)
