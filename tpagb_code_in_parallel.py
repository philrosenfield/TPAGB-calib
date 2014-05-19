from IPython import parallel
import time
import matplotlib.pyplot as plt
import numpy as np
import sfh_tests_multi_proc as sfh_tests
import ResolvedStellarPops as rsp
import stats
import model_plots
from data_plots import plot_cum_sum_sfr, plot_cmd_lf

def caller(vSFH, vsfh_kws):
    return vSFH.vary_the_SFH(**vsfh_kws)

def main(targets, cmd_inputs, nsfhs, dry_run=False, comp_corr=False):
    '''
    calls sfh_tests_multi_proc.sfh_tests_multi_proc in most basic way possible
    for up to 2 * available processors. Target & cmd_inputs are distributed and
    nsfhs are all done per processor. So with 12 processors, you can do up to
    24 target and cmd_inputs in one call, or code something better for a change.
    '''
    clients = parallel.Client()
    clients.block = False

    clients[:].execute('cd ~/research/TP-AGBcalib/code/TPAGB-calib/')
    clients[:].execute('import sfh_tests_multi_proc as sfh_tests')

    if comp_corr is True:
        table_file = 'comp_corr'
    else:
        table_file = 'default'
    vSFHs, vsfh_kws = sfh_tests.prepare_vsfh_run(targets, cmd_inputs, nsfhs,
                                                 dry_run=dry_run,
                                                 table_file=table_file)

    if dry_run is True:
        timeout = 10
    else:
        timeout = 900

    # find a better way to run this all at once, what if I need three times
    # through?
    nprocs = len(clients)
    nvsfhs = len(vSFHs)
    ntimes = np.min([nprocs, nvsfhs])
    ndiff = np.abs(nvsfhs - nprocs)

    if ndiff > nprocs:
        print 'need a for loop, too many processes code code code man'
        import sys
        sys.exit()

    print 'calling first set'
    res = [clients[i].apply(caller, vSFHs[i], vsfh_kws,)
           for i in range(ntimes)]

    while False in [r.ready() for r in res]:
        print 'waiting on first set'
        time.sleep(timeout)
        print 'checking first set...'

    print 'writing first set'
    [vSFHs[i].write_results(res[i].result) for i in range(ntimes)]

    print 'calling second set'
    res2 = [clients[i].apply(caller, vSFHs[i+ntimes], vsfh_kws,)
            for i in range(ndiff)]

    while False in [r.ready() for r in res2]:
        print 'waiting on second set'
        time.sleep(timeout)
        print 'checking second set...'

    print 'writing second set'
    [vSFHs[i+ntimes].write_results(res2[i].result) for i in range(ndiff)]

    print 'done.'

def lf_figs(targets, cmd_inputs, nsfhs, outfile_dir='default', extra_str='',
            default_kw=None, comp_corr=False, example=True):
    import os
    if comp_corr is True:
        table_file = 'comp_corr'
    else:
        table_file = 'default'
    vSFHs, vsfh_kws = sfh_tests.prepare_vsfh_run(targets, cmd_inputs, nsfhs,
                                                 vsfh_kw={'outfile_loc': outfile_dir,
                                                          'extra_str': extra_str,
                                                          'table_file': table_file},
                                                 default_kw=default_kw)

    if comp_corr is True:
        extra_str += '_comp'

    for i in range(len(vSFHs)):
        pl = sfh_tests.Plotting(vSFHs[i])
        pl.compare_to_gal(extra_str=extra_str,
                          completeness_correction=comp_corr)
        # example LF from the model
    if example is True:
        for i in range(len(vSFHs)):
            pl = sfh_tests.Plotting(vSFHs[i])
            best = rsp.fileIO.get_files(os.path.join(outfile_dir, vSFHs[i].target,
                                                      vSFHs[i].agb_mod, 'mc'),
                                         '*best.dat')
            if len(best) == 0:
                continue
            pl.compare_to_gal(narratio=False, add_stage_lfs='all',
                              extra_str='no_data', plot_data=False,
                              completeness_correction=comp_corr,
                              plot_models=False,
                              trilegal_output=best[0])

    return

def narratio_table(outfile_dir):
    narratio_files = rsp.fileIO.get_files(outfile_dir, '*narratio*dat')
    stats.narratio_table(narratio_files)
    return

def chi2_stats(targets, cmd_inputs, outfile_dir='default', extra_str=''):
    chi2_files = stats.write_chi2_table(targets, cmd_inputs,
                                        outfile_loc=outfile_dir,
                                        extra_str=extra_str)
    chi2_dicts = stats.result2dict(chi2_files)
    stats.chi2plot(chi2_dicts, outfile_loc=outfile_dir)
    chi2_files = stats.write_chi2_table(targets, cmd_inputs,
                                        outfile_loc=outfile_dir,
                                        extra_str=extra_str,
                                        just_gauss=True)
    plt.close('all')
    return


def analysis(targets, cmd_inputs, nsfhs, outfile_dirs, extra_str='',
             comp_corr=False, default_kw=None):
    default_kw = default_kw or {}
    # _mstar or _cstar:
    if 'star' in extra_str:
        from TPAGBparams import snap_src
        default_kw = dict({'galaxy_input_src': snap_src + '/input/tests/',
                           'galaxy_input_search_fmt': '*%s' + '%s.dat' % extra_str}.items()
                          + default_kw.items())

    for out_dir in outfile_dirs:
        lf_figs(targets, cmd_inputs, nsfhs, outfile_dir=out_dir,
                extra_str=extra_str, default_kw=default_kw,
                comp_corr=comp_corr)
        #chi2_stats(targets, cmd_inputs, outfile_dir=out_dir,
        #           extra_str=extra_str)
        #narratio_table(out_dir)
    #model_plots.plot_random_sfhs(targets)

def all_data_plots(targets):
    plot_cum_sum_sfr(targets)
    [[plot_cmd_lf(target, band) for target in targets] for band in ['opt', 'ir']]

if __name__ == '__main__':
    # could be an input file:
    targets = ['scl-de1', 'ddo71', 'kkh37', 'ddo78', 'hs117', 'ngc2976-deep']
    #targets = ['scl-de1']
    cmd_inputs = ['cmd_input_CAF09_S_NOV13.dat',
                  'cmd_input_CAF09_S_NOV13eta0.dat',
                  'cmd_input_CAF09_S_OCT13.dat']
    models = [c.split('S_')[1].replace('.dat', '') for c in cmd_inputs]
    #outfile_dirs = ['/home/rosenfield/research/TP-AGBcalib/SNAP/models/varysfh/match-hmc/']
    outfile_dirs = ['/home/rosenfield/research/TP-AGBcalib/SNAP/models/varysfh/comp_corr']
    nsfhs = 5
    dry_run = True
    comp_corr = True

    #main(targets, cmd_inputs, nsfhs, dry_run=dry_run, comp_corr=comp_corr)
    #all_data_plots(targets)
    #model_plots.agb_lifetimes(models, z='all')
    #model_plots.compare_agb_lifetimes()
    analysis(targets, cmd_inputs, nsfhs, outfile_dirs, comp_corr=comp_corr)
    #model_plots.tpagb_mass_histograms(chi2_location=outfile_dirs[0],
    #                                   band='opt', dry_run=True, model='nov13',
    #                                   model_src=outfile_dirs[0], force=True,
    #                                   cumsum=False)
    #[model_plots.compare_mass_loss(masses=m, z=0.001, paola=True)
    # for m in [1., 2.]]
