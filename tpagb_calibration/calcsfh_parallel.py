"""Run calcsfh or hybridMC in Parallel (using subprocess)"""
import argparse
import logging
import os
import subprocess
import sys

import numpy as np

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# Could be in a config or environ
calcsfh = '$HOME/research/match2.5/bin/calcsfh'
zcombine = '$HOME/research/match2.5/bin/zcombine'
hybridmc = '$HOME/research/match2.5/bin/hybridMC'

def test_files(prefs, run_calcsfh=True):
    """make sure match input files exist"""
    return_code = 0
    for pref in prefs:
        if run_calcsfh:
            pfiles = calcsfh_existing_files(pref)
        else:
            pfiles = [hybridmc_existing_files(pref)]
        test = [os.path.isfile(f) for f in pfiles]
        if False in test:
            logger.error('missing a file in {}'.format(pref))
            logger.error(pfiles)
            return_code += 1
    if return_code > 0:
        sys.exit(2)
    return


def calcsfh_existing_files(pref):
    """file formats for param match and matchfake"""
    pref = pref.strip()
    param = pref + '.param'
    match = pref + '.match'
    fake = pref + '.matchfake'
    return (param, match, fake)


def calcsfh_new_files(pref):
    """file formats for match grid, sdout, and sfh file"""
    pref = pref.strip()
    out = pref + '.out'
    scrn = pref + '.scrn'
    sfh = pref + '.sfh'
    return (out, scrn, sfh)


def hybridmc_existing_files(pref):
    """file formats for the HMC, based off of calcsfh_new_files"""
    pref = pref.strip()
    mcin = pref + '.out.dat'
    return mcin


def hybridmc_new_files(pref):
    """file formats for HybridMC output and the following zcombine output"""
    pref = pref.strip()
    mcmc = pref + '.mcmc'
    mcscrn = mcmc + '.scrn'
    mczc = mcmc + '.zc'
    return (mcmc, mcscrn, mczc)


def run_parallel(prefs, dry_run=False, nproc=8, run_calcsfh=True):
    """run calcsfh and zcombine in parallel, flags are hardcoded."""
    test_files(prefs, run_calcsfh)

    # calcsfh
    cmd1 = '{0} {1} {2} {3} {4} -PARSEC -mcdata -kroupa -zinc -sub=v2 > {5}'
    # zcombine
    cmd2 = '{0} {1} -bestonly > {2}'
    # hybridmc
    cmd3 = '{0} {1} {2} -tint=2.0 -nmc=10000 -dt=0.015 > {3}'
    # zcombine w/ hybrid mc
    cmd4 = '{0} {1} -unweighted -medbest -jeffreys -best={2}'

    niters = np.ceil(len(prefs) / float(nproc))
    sets = np.arange(niters * nproc, dtype=int).reshape(niters, nproc)
    logging.debug('{} prefs, {} niters'.format(len(prefs), niters))

    for j, iset in enumerate(sets):
        # don't use not needed procs
        iset = iset[iset < len(prefs)]
        
        # run calcsfh
        procs = []
        for i in iset:
            if run_calcsfh:
                param, match, fake = calcsfh_existing_files(prefs[i])
                out, scrn, sfh = calcsfh_new_files(prefs[i])
                cmd = cmd1.format(calcsfh, param, match, fake, out, scrn)
            else:
                mcin = hybridmc_existing_files(prefs[i])
                mcmc, mcscrn, mczc = hybridmc_new_files(prefs[i])
                cmd = cmd3.format(hybridmc, mcin, mcmc, mcscrn)
            if not dry_run:
                procs.append(subprocess.Popen(cmd, shell=True))
            logger.info(cmd)
        
        # wait for calcsfh
        if not dry_run:
            [p.wait() for p in procs]
            logger.debug('calcsfh or hybridMC set {} complete'.format(j))
        
        # run zcombine
        procs = []
        for i in iset:
            if run_calcsfh:
                out, scrn, sfh = calcsfh_new_files(prefs[i])
                zcom = cmd2.format(zcombine, out, sfh)
            else:
                zcom = cmd4.format(zcombine, mcmc, mczc)
            if not dry_run:
                procs.append(subprocess.Popen(zcom, shell=True))
            logger.info(zcom)
        
        # wait for zcombine
        if not dry_run:
            [p.wait() for p in procs]
            logger.debug('zcombine set {} complete'.format(j))


def main(argv):
    """parse in put args, setup logger, and call run_parallel"""
    parser = argparse.ArgumentParser(description="Run calcsfh in parallel. Note: bg cmd, if in use, need to be in the current folder")

    parser.add_argument('-d', '--dry_run', action='store_true',
                        help='only print commands')

    parser.add_argument('-v', '--verbose', action='store_true',
                        help='set logging to debug')

    parser.add_argument('-n', '--nproc', type=int, default=8,
                        help='number of processors')

    parser.add_argument('-m', '--hmc',  action='store_false',
                        help='run hybridMC (must be after a calcsfh run)')

    parser.add_argument('-f', '--logfile', type=str, default='calcsfh_parallel.log',
                        help='log file name')

    parser.add_argument('pref_list', type=argparse.FileType('r'),
                        help="list of prefixs to run on. E.g., ls */*.match | sed 's/.match//' > pref_list")

    args = parser.parse_args(argv)
    prefs = args.pref_list.readlines()
    
    handler = logging.FileHandler(args.logfile)
    if args.verbose:
        handler.setLevel(logging.DEBUG)
    else:
        handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    logger.info('running on {}'.format(', '.join([p.strip() for p in prefs])))
    run_parallel(prefs, dry_run=args.dry_run, nproc=args.nproc, run_calcsfh=args.hmc)


if __name__ == '__main__':
    main(sys.argv[1:])