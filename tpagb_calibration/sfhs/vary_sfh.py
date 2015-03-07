"""
Run many trilegal simulations and cull scaled LF functions to compare with data
"""
import argparse
import logging
import numpy as np
import os
import sys
import time

import ResolvedStellarPops as rsp

from IPython import parallel
from .star_formation_histories import StarFormationHistories

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

__all__ = ['VarySFHs']

def initialize_inputs():
    return {'file_origin': None,
            'filter1': None,
            'filter2': None,
            'galaxy_input': None,
            'outfile_loc': os.getcwd(),
            'target': None}


def jobwait(line=''):
    line += "\nfor job in `jobs -p`\ndo\n    echo $job\n    wait $job\ndone\n\n"
    return line

class VarySFHs(StarFormationHistories):
    '''
    run several variations of the age sfr z from MATCH SFH to produce
    simulated CMDs, LFs, and nagb/nrgb ratios.

    because the trilegal output can be so huge, there are many constraints
    needed to cull the model output down to what is necessary for the
    analysis.
    '''
    def __init__(self, inp_obj=None, input_file=None):
        '''
        galaxy_input is a template.
        '''
        # load SFH instance to make lots of trilegal runs
        self.input_file = input_file
        if input_file is not None:
            indict = rsp.fileio.load_input(input_file)
            indict = dict(initialize_inputs.items() + indict.items())
        if inp_obj is not None:
            indict = inp_obj.__dict__

        if inp_obj.nsfhs > 1:
            StarFormationHistories.__init__(self, inp_obj.hmc_file,
                                            inp_obj.file_origin)

        [self.__setattr__(k, v) for k, v in indict.items()]

        cmd_input_file = os.path.split(self.cmd_input_file)[1]
        self.agb_mod = \
            cmd_input_file.replace('.dat', '').lower().replace('cmd_input_', '')

    def prepare_trilegal_output(self):
        # trilegal output format
        self.tname = os.path.join(self.outfile_loc,
                             'output_%s_%s_%s_%s' % (self.target, self.filter1,
                                                     self.filter2, self.agb_mod))
        self.triout_fmt = self.tname + '_%003i.dat'
        
    def prepare_galaxy_input(self, object_mass=None, overwrite=False):
        '''
        write the galaxy input file from a previously written template.
        simply overwrites the filename line to link to the new sfr
        file.
        '''
        self.galaxy_inputs = []
        galaxy_input = self.galaxy_input
        ext = '.' + galaxy_input.split('.')[-1]
        lines = open(galaxy_input).readlines()
        # line that links to sfr file.
        extra = ' '.join(lines[-3].split(' ')[1:])

        if object_mass is not None:
            extra2 = ' '.join(lines[-6].split()[1:])
        
        for i in range(len(self.sfr_files)):
            lines[-3] = ' '.join([self.sfr_files[i], extra])
            if object_mass is not None:
                lines[-6] = ' '.join(['%.4e' % object_mass, extra2]) + '\n'
            new_name = os.path.split(galaxy_input)[1].replace(ext, '_%003i' % i + ext)
            new_out = os.path.join(self.outfile_loc, new_name)
            if not os.path.isfile(new_out) or overwrite:
                with open(new_out, 'w') as f:
                    f.write(''.join(lines))
                logger.info('wrote {}'.format(new_out))
            else:
                logger.info('not overwritting {}'.format(new_out))
            self.galaxy_inputs.append(new_out)

    def vary_the_SFH(self, random_sfr=True, random_z=False,
                     zdisp=False, overwrite=False, object_mass=None):
        '''make the sfhs, make the galaxy inputs'''
        new_fmt = '{}{}_tri_%003i.sfr'.format(self.target, self.filter1) 
        outfile_fmt = os.path.join(self.outfile_loc, new_fmt)
        self.sfr_files = self.make_many_trilegal_sfhs(nsfhs=self.nsfhs,
                                                      outfile_fmt=outfile_fmt,
                                                      random_sfr=random_sfr,
                                                      random_z=random_z,
                                                      zdisp=zdisp,
                                                      dry_run=overwrite)

        self.prepare_galaxy_input(overwrite=overwrite, object_mass=object_mass)

        return

    def run_once(self, galaxy_input=None, triout=None, ite=0, overwrite=False):
        """call trilegal and convert the output file to hdf5"""
        import ResolvedStellarPops as rsp
        flag = 0
        ver = 2.3
        call = ''
        #print('cmd: {} galinp: {} out: {} dryrun: {}'.format(self.cmd_input_file,
        #                                                     galaxy_input,
        #                                                     triout,
        #                                                     dry_run))
        if os.path.isfile(triout) and not overwrite:
            logger.warning('{} exists, will overwrite if no hdf5 file found'.format(triout))
            #flag += 1

        hdf5file = rsp.fileio.replace_ext(triout, 'hdf5')
        if os.path.isfile(hdf5file) and not overwrite:
            logger.warning('{} already exists, not calling trilegal'.format(hdf5file))
            flag += 1
        
        if flag < 1:
            call = 'nice -n +19 taskset -c {0} code_{1}/main'.format(ite, ver)
            call += ' -f {0} -a -l {1} {2} > {2}.scrn'.format(self.cmd_input_file,
                                                               galaxy_input,
                                                               triout)
        return call

    def call_run(self, dry_run=False, nproc=8, overwrite=False):
        """Call run_once or run_parallel depending on self.nsfh value"""

        self.prepare_trilegal_output()

        if self.nsfhs <= 1:
            # don't run parallel.
            cmd = self.run_once(galaxy_input=self.galaxy_input,
                          triout=self.tname + '_bestsfr.dat',
                          overwrite=overwrite)
        else:
            cmd = self.run_many(nproc=nproc, overwrite=overwrite)
        return cmd

    def run_many(self, nproc=8, overwrite=False):
        """
        Call self.run_once a bunch of times
        """
        self.vary_the_SFH(random_sfr=True, random_z=False, zdisp=False,
                          overwrite=overwrite, object_mass=None)

        # How many sets of calls to the max number of processors
        niters = np.ceil(self.nsfhs / float(nproc))
        sets = np.arange(niters * nproc, dtype=int).reshape(niters, nproc)

        line = ''
        for j, iset in enumerate(sets):
            # don't use not needed procs
            iset = iset[iset < self.nsfhs]
            for i in range(len(iset)):
                cmd = self.run_once(galaxy_input=self.galaxy_inputs[iset[i]],
                                    triout=self.triout_fmt % iset[i], ite=i,
                                    overwrite=overwrite)
                line += '{}\n'.format(cmd)
            line += jobwait()

        return line

def call_VarySFH(input_file, loud=False, nproc=8, outfile=None,
                 overwrite=False):
    # set up logging
    from IPython.config import Application
    logger = Application.instance().log
    handler = logging.FileHandler('{}_vary_sfh.log'.format(input_file))
    logger.setLevel(logging.DEBUG)
    if loud:
        handler.setLevel(logging.DEBUG)
    else:
        handler.setLevel(logging.DEBUG)
    formatter = logging.Formatter('%(asctime)s - %(name)s - %(levelname)s - %(message)s')
    handler.setFormatter(formatter)
    logger.addHandler(handler)
    logger.info('logger writing to {}_vary_sfh.log'.format(input_file))
    # set up input parameters
    inp_obj = rsp.fileio.InputParameters(default_dict=initialize_inputs())
    inp_obj.input_file = input_file
    inp_obj.add_params(rsp.fileio.load_input(inp_obj.input_file), loud=loud)

    #  do it!
    vsh = VarySFHs(inp_obj=inp_obj)
    #import pdb; pdb.set_trace()
    line = vsh.call_run(nproc=nproc, overwrite=overwrite)
    if outfile is None:
        print(line)
    else:
        logger.info('output file: {}'.format(outfile))
        with open(outfile, 'a') as out:
            out.write(line)
    return

def main(argv):
    """
    call vsfh.run_parallel with command line options and set up logger.
    """
    parser = argparse.ArgumentParser(description="Run trilegal many times by \
                                     randomly sampling SFH uncertainies")

    parser.add_argument('-v', '--verbose', action='store_true',
                        help='verbose mode')

    parser.add_argument('-n', '--nproc', type=int, default=8,
                        help='number of processors')

    parser.add_argument('-o', '--outfile', type=str, default='trilegal_script.sh',
                        help='name of output script')

    parser.add_argument('-f', '--overwrite', action='store_true',
                        help='write call to trilegal even if output file exists')

    parser.add_argument('name', type=str, help='input file')

    args = parser.parse_args(argv)

    call_VarySFH(args.name, loud=args.verbose, nproc=args.nproc,
                 outfile=args.outfile, overwrite=args.overwrite)


if __name__ == '__main__':
    main(sys.argv[1:])
