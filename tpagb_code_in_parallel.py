from IPython import parallel
import sfh_tests_multi_proc as sfh_tests
import time
import numpy as np

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

if __name__ == '__main__':
    # could be an input file:
    targets = ['scl-de1', 'ddo71', 'kkh37', 'ddo78', 'hs117', 'ngc2976-deep']
    cmd_inputs = ['cmd_input_CAF09_S_NOV13.dat',
                  'cmd_input_CAF09_S_NOV13eta0.dat',
                  'cmd_input_CAF09_S_OCT13.dat']
    nsfhs = 50
    dry_run = False
    comp_corr = True

    main(targets, cmd_inputs, nsfhs, dry_run=dry_run, comp_corr=comp_corr)