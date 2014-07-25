from IPython import parallel
import time
import numpy as np


def caller(vsfh, vsfh_kw={}):
    return vsfh.vary_the_SFH(**vsfh_kw)

def main(vsfhs, vsfh_kws=None, timeout=900, dry_run=False):
    '''
    calls sfh_tests_multi_proc.sfh_tests_multi_proc in most basic way possible
    for up to 2 * available processors. Target & cmd_inputs are distributed and
    nsfhs are all done per processor. So with 12 processors, you can do up to
    24 target and cmd_inputs in one call, or code something better for a change.
    '''
    clients = parallel.Client()
    clients.block = False

    #clients[:].execute('cd ~/research/TP-AGBcalib/code/TPAGB-calib/tpagb_calibration')
    #clients[:].execute('from sfhs.vary_sfh import ')

    if dry_run is True:
        timeout = 10

    # find a better way to run this all at once, what if I need three times
    # through?
    nprocs = len(clients)
    nvsfhs = len(vsfhs)
    ntimes = np.min([nprocs, nvsfhs])
    ndiff = np.abs(nvsfhs - nprocs)

    if vsfh_kws is None:
        vsfh_kws = [{}] * nvsfhs

    if ndiff > nprocs:
        print 'need a for loop, too many processes code code code man'
        import sys
        sys.exit()

    print 'calling first set'
    res = [clients[i].apply(caller, vsfhs[i], vsfh_kws[i],)
           for i in range(ntimes)]

    while False in [r.ready() for r in res]:
        print 'waiting on first set'
        time.sleep(timeout)
        print 'checking first set...'

    #print 'writing first set'
    #[vsfhs[i].write_results(res[i].result) for i in range(ntimes)]

    print 'calling second set'
    res2 = [clients[i].apply(caller, vsfhs[i+ntimes], vsfh_kws[i+ntimes],)
            for i in range(ndiff)]

    while False in [r.ready() for r in res2]:
        print 'waiting on second set'
        time.sleep(timeout)
        print 'checking second set...'

    #print 'writing second set'
    #[vsfhs[i+ntimes].write_results(res2[i].result) for i in range(ndiff)]

    print 'done.'
