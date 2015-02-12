import argparse
import logging
import os
import subprocess
import sys

import numpy as np

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

calcsfh = '$HOME/research/match2.5/bin/calcsfh'
zcombine = '$HOME/research/match2.5/bin/zcombine'

def test_files(prefs):
    """make sure match input files exist"""
    return_code = 0
    for pref in prefs:
        p, m, f = existing_files(pref)
        test = [os.path.isfile(f) for f in [p, m, f]]
        if False in test:
            logger.error('missing a file in {}'.format(pref))
            logger.error('{} {} {}'.format(p, m, f))
            return_code += 1
    if return_code > 0:
        sys.exit(2)
    return


def existing_files(pref):
    """file formats for param match and matchfake"""
    pref = pref.strip()
    param = pref + '.param'
    match = pref + '.match'
    fake = pref + '.matchfake'
    return (param, match, fake)


def new_files(pref):
    """file formats for match grid, sdout, and sfh file"""
    pref = pref.strip()
    out = pref + '.out'
    scrn = pref + '.scrn'
    sfh = pref + '.sfh'
    return (out, scrn, sfh)


def run_parallel(prefs, dry_run=False, nproc=8):
    """
    run calcsfh and zcombine in parallel, flags are currently hardcoded.
    """
    test_files(prefs)
    
    # calcsfh
    cmd1 = '{0} {1} {2} {3} {4} -kroupa -zinc > {5}'
    # zcombine
    cmd2 = '{0} {1} -bestonly > {2}'

    niters = np.ceil(len(prefs) / float(nproc))
    sets = np.arange(niters * nproc, dtype=int).reshape(niters, nproc)
    logging.debug('{} prefs, {} niters'.format(len(prefs), niters))

    for j, iset in enumerate(sets):
        # don't use not needed procs
        iset = iset[iset < len(prefs)]
        
        # run calcsfh
        procs = []
        for i in iset:
            param, match, fake = existing_files(prefs[i])
            out, scrn, sfh = new_files(prefs[i])
            csfh = cmd1.format(calcsfh, param, match, fake, out, scrn)
            procs.append(subprocess.Popen(csfh, shell=True))
            logger.debug(csfh)
        
        # wait for calcsfh
        [p.wait() for p in procs]
        
        # run zcombine
        procs = []
        for i in iset:
            out, scrn, sfh = new_files(prefs[i])
            zcom = cmd2.format(zcombine, out, sfh)
            procs.append(subprocess.Popen(zcom, shell=True))
            logger.debug(zcom)
        
        # wait for zcombine
        [p.wait() for p in procs]
        
        logger.debug('set {} complete'.format(j))


def main(argv):
    parser = argparse.ArgumentParser(description="Run calcsfh in parallel")

    parser.add_argument('-d', '--dry_run', action='store_true',
                        help='only print commands')

    parser.add_argument('-v', '--verbose', action='store_true',
                        help='set logging to debug')

    parser.add_argument('-n', '--nproc', type=int, default=8,
                        help='number of processors')

    parser.add_argument('pref_list', type=argparse.FileType('r'),
                        help="list of prefixs, to make try: ls */*.match | sed 's/.match//' > pref_list recal that bg cmds, if in use, need to be in the current folder")

    args = parser.parse_args(argv)
    prefs = args.pref_list.readlines()
    
    handler = logging.FileHandler('calcsfh_parallel.log')
    if args.verbose:
        handler.setLevel(logging.DEBUG)
    else:
        handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)

    run_parallel(prefs, dry_run=args.dry_run, nproc=args.nproc)


if __name__ == '__main__':
    main(sys.argv[1:])