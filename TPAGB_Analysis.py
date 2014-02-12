import sfh_tests_multi_proc
import ResolvedStellarPops as rsp
import stats
from TPAGBparams import snap_src
import model_plots

# nsfhs dont matter at all here but need to be defined for vsfh...
def full_run():
    targets = ['ddo71', 'hs117', 'kkh37', 'ngc2976-deep', 'scl-de1', 'ddo78']
    cmd_inputs = ['cmd_input_CAF09_S_NOV13.dat',
                  'cmd_input_CAF09_S_NOV13eta0.dat',
                  'cmd_input_CAF09_S_OCT13.dat']
    nsfhs = 50
    outfile_dirs = ['/home/rosenfield/research/TP-AGBcalib/SNAP/models/varysfh/02092014/']
    extra_str = ''
    return targets, cmd_inputs, nsfhs, outfile_dirs, extra_str

def feb_run():
    targets = ['hs117', 'ngc2976-deep']
    cmd_inputs = ['cmd_input_CAF09_S_FEB14.dat']
    nsfhs = 50
    outfile_dirs = ['/home/rosenfield/research/TP-AGBcalib/SNAP/models/varysfh/feb_test/']
    extra_str = ''
    return targets, cmd_inputs, nsfhs, outfile_dirs, extra_str
    
def cstar():
    targets = ['ngc2976-deep', 'scl-de1']
    cmd_inputs = ['cmd_input_CAF09_S_NOV13.dat',
                  'cmd_input_CAF09_S_NOV13eta0.dat',
                  'cmd_input_CAF09_S_OCT13.dat']
    nsfhs = 50
    outfile_dirs = ['/home/rosenfield/research/TP-AGBcalib/SNAP/models/cstar_only']
    extra_str = '_cstar'
    return targets, cmd_inputs, nsfhs, outfile_dirs, extra_str

def mstar():
    targets = ['ngc2976-deep', 'scl-de1']
    cmd_inputs = ['cmd_input_CAF09_S_NOV13.dat',
                  'cmd_input_CAF09_S_NOV13eta0.dat',
                  'cmd_input_CAF09_S_OCT13.dat']
    nsfhs = 50
    outfile_dirs = ['/home/rosenfield/research/TP-AGBcalib/SNAP/models/mstar_only']
    extra_str = '_mstar'
    return targets, cmd_inputs, nsfhs, outfile_dirs, extra_str

targets, cmd_inputs, nsfhs, outfile_dirs, extra_str = mstar()
default_kw = {'galaxy_input_src': snap_src + '/input/tests/',
              'galaxy_input_search_fmt': '*%s' + '%s.dat' % extra_str}

targets, cmd_inputs, nsfhs, outfile_dirs, extra_str = feb_run()
default_kw = None

vsfh_kw = {'extra_str' : extra_str}

def all_stats(targets, cmd_inputs, nsfhs, outfile_dir='default', extra_str='', default_kw=None):
    vSFHs, vsfh_kws = sfh_tests_multi_proc.prepare_vsfh_run(targets, cmd_inputs, nsfhs,
                                                            vsfh_kw={'outfile_loc': outfile_dir,
                                                                     'extra_str': extra_str},
                                                            default_kw=default_kw)
    chi2_files = stats.write_chi2_table(targets, cmd_inputs, outfile_loc=outfile_dir, extra_str=extra_str)
    narratio_files = rsp.fileIO.get_files(outfile_dir, '*narratio*dat')
    chi2_dicts = stats.result2dict(chi2_files)
    stats.narratio_table(narratio_files)
    #axs = stats.chi2plot(chi2_dicts, outfile_loc=outfile_dir)
    axs = stats.chi2plot2(chi2_dicts, outfile_loc=outfile_dir)
    chi2_files = stats.write_chi2_table(targets, cmd_inputs, outfile_loc=outfile_dir, extra_str=extra_str,
                                        just_gauss=True)
    pl = [sfh_tests_multi_proc.Plotting(v) for v in vSFHs]
    _ = [pl[i].compare_to_gal(extra_str=extra_str) for i in range(len(pl))]
    #plt.close('all')
    return


_ = [all_stats(targets, cmd_inputs, nsfhs, outfile_dir=out_dir, extra_str=extra_str,
               default_kw=default_kw) for out_dir in outfile_dirs]

#model_plots.compare_mass_loss(mass=1, z=0.001)

#plt.close('all')

#model_plots.agb_lifetimes(['NOV13'], z='all')


