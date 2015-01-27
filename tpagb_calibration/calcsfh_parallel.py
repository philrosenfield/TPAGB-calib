import argparse
import logging
import os
import shlex
import subprocess
import sys
import time

from IPython.config import Application
from IPython import parallel
import numpy as np


logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

calcsfh = '$HOME/research/match2.5/bin/calcsfh'
zcombine = '$HOME/research/match2.5/bin/zcombine'

def test_files(prefs):
    """
    make sure match input files exist
    """
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
    pref = pref.strip()
    param = pref + '.param'
    match = pref + '.match'
    fake = pref + '.matchfake'
    return (param, match, fake)
    
def new_files(pref):
    pref = pref.strip()
    out = pref + '.out'
    scrn = pref + '.scrn'
    sfh = pref + '.sfh'
    return (out, scrn, sfh)


def run_once(pref, dry_run=False):
    param, match, fake = existing_files(pref)
    out, scrn, sfh = new_files(pref)
    cmd1 = '{0} {1} {2} {3} {4} -kroupa -zinc > {5}'.format(calcsfh, param,
                                                            match, fake, out,
                                                            scrn)
    cmd2 = '{0} {1} -bestonly > {2}'.format(zcombine, out, sfh)
    logger.info(cmd1)
    logger.info(cmd2)

    #if not dry_run:
    p = subprocess.Popen(cmd1, shell=True)
    p.wait()
    q = subprocess.Popen(cmd2, shell=True)
    q.wait()
    return


def run_parallel(prefs, dry_run=False, nproc=8, start=45):
    """
    """
    def setup_parallel():
        """
        I would love a better way to do this.
        """
        clients = parallel.Client()
        clients.block = False
        clients[:].use_dill()
        clients[:].execute('import os')
        clients[:].execute('import logging')
        clients[:].execute('import numpy as np')
        clients[:].execute('from IPython.config import Application')
        clients[:]['run_once'] = run_once
        clients[:]['new_files'] = new_files
        clients[:]['existing_files'] = existing_files
        clients[:].execute('logger = Application.instance().log')
        return clients

    test_files(prefs)

    if len(prefs) == 1:
        run_once(prefs[0], dry_run=dry_run)
        return

    try:
        clients = parallel.Client()
    except IOError:
        logger.debug('Starting ipcluster. Waiting {}s for spin up'.format(start))
        os.system('ipcluster start --n={} &'.format(nproc))
        time.sleep(start)

    # find looping parameters. How many sets of calls to the max number of
    # processors

    niters = np.ceil(len(prefs) / float(nproc))
    sets = np.arange(niters * nproc, dtype=int).reshape(niters, nproc)
    logging.debug('{} prefs, {} niters, {} sets'.format(len(prefs), niters, sets))
    
    # in case it takes more than start sec to spin up clusters, set up as
    # late as possible
    clients = setup_parallel()

    for j, iset in enumerate(sets):
        # don't use not needed procs
        iset = iset[iset < len(prefs)]

        # parallel call to run
        res = [clients[i].apply(run_once, prefs[i], dry_run,)
               for i in range(len(iset))]

        logger.debug('waiting on set {} of {}'.format(j, niters))
        while False in [r.ready() for r in res]:
            time.sleep(1)
        logger.debug('set {} complete'.format(j))
    os.system('ipcluster stop')

def main(argv):
    parser = argparse.ArgumentParser(description="Run calcsfh in parallel")

    parser.add_argument('-d', '--dry_run', action='store_true',
                        help='only print commands')

    parser.add_argument('-v', '--verbose', action='store_true',
                        help='set logging to debug')

    parser.add_argument('-n', '--nproc', type=int, default=8,
                        help='number of processors')

    parser.add_argument('pref_list', type=argparse.FileType('r'),
                        help="list of prefixs, try: ls */*match | sed 's/.match//' > pref_list")

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