"""
Add AST corrections to trilegal catalogs
Note: single catalog mode is untested
"""
import argparse
import logging
import os
import sys

import ResolvedStellarPops as rsp

from ..TPAGBparams import snap_src

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

# where the matchfake files live
matchfake_loc = os.path.join(snap_src, 'data', 'galaxies')


def load_trilegal_catalog(trilegal_catalog):
    '''read a table into a SimGalaxy object'''
    return rsp.SimGalaxy(trilegal_catalog)


def make_ast_corrections(trilegal_catalogs, target, outfiles='default',
                         diag_plot=False, overwrite=True):
    """
    apply ast corrections from fake files found in matchfake_loc/*[target]*
    see rsp.ast_correct_starpop
    """
    if type(outfiles) is str:
        outfmt = 'default'
    else:
        outfmt = 'supplied'

    # search string for fake files
    search_str = '*{}*.matchfake'.format(target.upper())
    
    fakes = rsp.fileio.get_files(matchfake_loc, search_str)
    logger.info('fake files found: {}'.format(fakes))
    asts = [rsp.ASTs(f) for f in fakes]
    
    for i, trilegal_catalog in enumerate(trilegal_catalogs):
        logger.info('working on {}'.format(trilegal_catalog))
        sgal = load_trilegal_catalog(trilegal_catalog)
        # "overwrite" (append columns) to the existing catalog by default
        if outfmt == 'default':
            outfile = trilegal_catalog
        else:
            outfile = outfiles[i]

        # do the ast corrections
        [rsp.ast_correct_starpop(sgal, asts_obj=ast, overwrite=overwrite,
                                 outfile=outfile, diag_plot=diag_plot)
         for ast in asts]
    return

def main(argv):
    """
    Make AST corrections to trilegal catalog(s)
    
    usage:
    python add_asts.py -vda ~/research/TP-AGBcalib/SNAP/varysfh/kkh37
    
    if the target directory name is different than it is in the matchfake file name:
    python analyze.py -vda -t ugc4459 ~/research/TP-AGBcalib/SNAP/varysfh/ugc-04459
    """
    parser = argparse.ArgumentParser(description="Cull useful info from \
                                                  trilegal catalog")

    parser.add_argument('-d', '--directory', action='store_true',
                        help='opperate on *.dat files in a directory')

    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose mode')

    parser.add_argument('-t', '--target', type=str, help='target name')

    parser.add_argument('name', type=str, nargs='*',
                        help='trilegal catalog or directory if -d flag')

    args = parser.parse_args(argv)


    if not args.target:   
        if args.directory:
            target = os.path.split(args.name[0])[1]
        else:
            target = tricat.split('_')[1]
    else:
        target = args.target

    # set up logging
    handler = logging.FileHandler('{}_analyze.log'.format(target))
    if args.verbose:
        handler.setLevel(logging.DEBUG)
    else:
        handler.setLevel(logging.INFO)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.info('using matchfake location: {}'.format(matchfake_loc))
    
    # assume trilegal was run with outfile extension == .dat
    if args.directory:
        tricats = rsp.fileio.get_files(args.name[0], '*dat')
    else:
        tricats = args.name
        
    if args.verbose:
        logger.info('working on target: {}'.format(target))
        
    make_ast_corrections(tricats, target)

    return

if __name__ == "__main__":
    main(sys.argv[1:])