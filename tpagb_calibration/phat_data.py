import itertools
import numpy as np
import os
import sys
import time

import calibrate
from plotting import Plotting
from pop_synth.stellar_pops import rgb_agb_regions
import ResolvedStellarPops as rsp
from sfhs.vary_sfh import VarySFHs

def load_data_files(inputs):
    data_src = os.path.join(inputs.base_dir, inputs.data_src)
    if type(inputs.opt_data) is str:
        inputs.opt_data = [inputs.opt_data]
    if type(inputs.ir_data) is str:
        inputs.ir_data = [inputs.ir_data]
    opt_files = [os.path.join(data_src, i) for i in inputs.opt_data]
    ir_files = [os.path.join(data_src, i) for i in inputs.ir_data]
    return opt_files, ir_files


def prepare_data_table(inputs):
    if hasattr(inputs, 'Mtrgb'):
        mtrgb = inputs.Mtrgb
    else:
        mtrgb = np.nan
    if not hasattr(inputs, 'ir_trgb'):
        inputs.ir_trgb = rsp.astronomy_utils.Mag2mag(inputs.Mtrgb, inputs.ir_filter2,
                                              inputs.photsys, dmod=inputs.dmod,
                                              Av=inputs.Av)
    if not hasattr(inputs, 'opt_trgb'):
        inputs.opt_trgb = rsp.astronomy_utils.Mag2mag(inputs.Mtrgb, inputs.opt_filter2,
                                               inputs.photsys, dmod=inputs.dmod,
                                               Av=inputs.Av)
    fmt = '%s %i %i %i %i %.2f %.2f %.2f %.2f \n'
    header = '# opt: F814W ir: F160W Mtrgb: %.2f dmod: %.2f Av: %.1f \n'
    header += '# opt, ir: color_min (%.1f, %.1f) trgb (%.2f, %.2f) '
    header += 'offsets (%.1f, %.1f) trgb exclude (%.1f, %.1f) \n'
    header += '# target nopt_rgb nopt_agb nir_rgb nir_agb opt_max ir_max'
    header += ' opt_min ir_min \n'

    out = open('phat_b21_data.dat', 'w')
    out.write(header % (mtrgb, inputs.dmod, inputs.Av,
                        inputs.opt_color_min, inputs.ir_color_min,
                        inputs.opt_trgb, inputs.ir_trgb, inputs.offsets[0], inputs.offsets[1],
                        inputs.trgb_excludes[0], inputs.trgb_excludes[1]))

    opt_files, ir_files = load_data_files(inputs)
    for i in range(len(opt_files)):
        opt_gal = rsp.Galaxy(opt_files[i], filetype=inputs.data_ftype, angst=False,
                             hla=False, filter1=inputs.opt_filter1,
                             filter2=inputs.opt_filter2)
        ir_gal = rsp.Galaxy(ir_files[i], filetype=inputs.data_ftype, angst=False,
                            hla=False, filter1=inputs.ir_filter1,
                            filter2=inputs.ir_filter2)
        opt_color_cut, = np.nonzero((opt_gal.color) > inputs.opt_color_min)
        ir_color_cut, = np.nonzero((ir_gal.color) > inputs.ir_color_min)
        ir_mag = ir_gal.mag2[ir_color_cut]
        opt_mag = opt_gal.mag2[opt_color_cut]
        opt_rgb, ir_rgb, opt_agb, ir_agb = \
            rgb_agb_regions(opt_gal, inputs.offsets, inputs.trgb_excludes,
                            inputs.opt_trgb, inputs.ir_trgb, opt_mag, ir_mag)
        out.write(fmt % (opt_files[i].split('.')[0], len(opt_rgb), len(opt_agb),
                         len(ir_rgb), len(ir_agb), opt_gal.mag2.max(),
                         ir_gal.mag2.max(), opt_gal.mag2.min(),
                         ir_gal.mag2.min()))


def load_numbers_table(table_file):
    dtype=[('target', '|S25'), ('nopt_rgb', '<f8'), ('nopt_agb', '<f8'),
           ('nir_rgb', '<f8'), ('nir_agb', '<f8'), ('opt_max', '<f8'),
           ('ir_max', '<f8'), ('opt_min', '<f8'), ('ir_min', '<f8')]
    return np.genfromtxt(table_file, dtype=dtype)


def prepare_vsfh_run(inputs, object_mass=None):
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
    cmd_input_files = inputs.cmd_input_files
    if type(cmd_input_files) is str:
        cmd_input_files = [cmd_input_files]

    targets = inputs.targets
    if type(targets) is str:
        targets = [targets]

    orig_outfile_loc = inputs.outfile_loc
    m31_data = load_numbers_table(os.path.join(inputs.base_dir, inputs.table_file))

    vsfh_kw = {'file_origin': inputs.sfh_file_origin,
               'nsfhs': inputs.nsfhs,
               'opt_filter1': inputs.opt_filter1,
               'opt_filter2': inputs.opt_filter2,
               'ir_filter1': inputs.ir_filter1,
               'ir_filter2': inputs.ir_filter2,
               'photsys': inputs.photsys,
               'opt_color_min': inputs.opt_color_min,
               'ir_color_min': inputs.ir_color_min,
               'dmod': inputs.dmod,
               'Av': inputs.Av,
               'Mtrgb': inputs.Mtrgb}
    vsfh_kws = []
    vsfhs = []

    for target, cmd_input in itertools.product(targets, cmd_input_files):
        agb_mod = cmd_input.replace('cmd_input_', '').lower().split('.')[0]

        galaxy_input = os.path.join(inputs.base_dir, inputs.galaxy_inp_src,
                                    target + inputs.galaxy_inp_ext)
        sfh_file = os.path.join(inputs.base_dir, inputs.sfh_src,
                                target + inputs.sfh_ext)

        if orig_outfile_loc == 'default':
            outfile_loc = os.path.join(inputs.base_dir, 'vary_sfh', target, agb_mod)

        rsp.fileio.ensure_dir(outfile_loc)

        data_line = m31_data['target'==target]
        opt_bins = np.arange(data_line['opt_min'], data_line['opt_max'], inputs.binsize)
        ir_bins = np.arange(data_line['ir_min'], data_line['ir_max'], inputs.binsize)

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


def phat_data(inputs, object_mass=None):
    vsfhs, vsfh_kws = prepare_vsfh_run(inputs)
    nruns = len(vsfhs) * inputs.nsfhs
    if nruns > 10:
        os.system('ipcluster start -n=%i' % inputs.nprocs)
        time.wait(45)
        calibrate.main(vsfhs, vsfh_kws=vsfh_kws)
        os.system('ipcluster stop')

    opt_files, ir_files = load_data_files(inputs)
    for i, vsfh in enumerate(vsfhs):
        if nruns < 10:
            if inputs.dry_run is False:
                res_dict = vsfh.vary_the_SFH(object_mass=object_mass, dry_run=inputs.dry_run)
                vsfh.write_results(res_dict)
        pl = Plotting(vsfh)
        gal_kw = {'filetype': inputs.data_ftype, 'angst': False, 'hla': False}
        opt_gal = rsp.Galaxy(opt_files[i], filter1='F475W', filter2='F814W',
                             **gal_kw)
        ir_gal = rsp.Galaxy(ir_files[i], filter1='F110W', filter2='F160W',
                            **gal_kw)

        pl.compare_to_gal(opt_gal, ir_gal, 28, 28, narratio=False, ylim=(.1, 1e5))


if __name__ == '__main__':
    import pdb; pdb.set_trace()
    inp_obj = rsp.fileio.InputFile(sys.argv[1])
    if hasattr(inp_obj, 'object_mass'):
        object_mass = inp_obj.object_mass
    else:
        object_mass = 1e7
    prepare_data_table(inp_obj)
    phat_data(inp_obj, object_mass=object_mass)