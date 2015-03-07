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

from .star_formation_histories import StarFormationHistories as SFH

logging.basicConfig(level=logging.DEBUG)
logger = logging.getLogger(__name__)

__all__ = ['VarySFHs']


def jobwait(line=''):
    """add bash script to wait for current jobs to finish"""
    line += "\nfor job in `jobs -p`\ndo\n    echo $job\n    wait $job\ndone\n\n"
    return line

class VarySFHs(SFH):
    '''run several variations of the age sfr z from SFH'''
    def __init__(self, inp_obj=None, input_file=None):
        """Vary the SFH from MATCH for a trilegal simulations
        
        Parameters
        ----------
            inp_obj : rsp.fileio.InputParameters object
                input parameters object
            input_file : path to file that can be read into a dictionary via rsp.fileio.load_input
        
            Necessary contents of input_file/inp_obj
            ------------------
            file_origin : str
                what type of SFH file (match-grid, match-hmc)
        
            filter1, filter2 : str, str
                V, I filters. Used only in file name conventions
        
            galaxy_input : str
                template galaxy input. object_mass and sfr-z file will be adjusted
        
            outfile_loc : str
                path to put the trilegal output files
        
            target : str
                name of observation target for file name conventions
        
            hmc_file : str
                path to the SFH file
        
            cmd_input_file : str
                path to the cmd input file to run TRILEGAL
        
            object_mass : str, will be converted to float
                optional, overwrite the mass set in galaxy_input
        
            nsfhs : str, will be converted to int
                number of sfhs to sample
        """
        # load SFH instance to make lots of trilegal runs
        self.input_file = input_file
        if input_file is not None:
            indict = rsp.fileio.load_input(input_file)

        if inp_obj is not None:
            indict = inp_obj.__dict__

        self.initialize_inputs()

        if self.nsfhs > 1:
            # load in hmc data to self.data
            SFH.__init__(self, self.hmc_file, self.file_origin)

        # setup file formats
        self.trilegal_file_fmt()

    def initialize_inputs(self):
        """load input parameters needed for vary_sfh"""
        # parameters needed
        inputs = ['file_origin', 'filter1', 'filter2', 'galaxy_input',
                  'outfile_loc', 'target', 'hmc_file', 'cmd_input_file',
                  'object_mass', 'nsfhs']

        needed = [k for k in inputs if not k in indict.keys()]
        if len(needed) > 0:
            logger.error('missing needed input parameters: {}'.format(needed))

        unused = [k for k in indict.keys() if not k in inputs]
        if len(unused) > 0:
            logger.warning('not using {}'.format(unused))

        [self.__setattr__(k, v) for k, v in indict.items() if k in inputs]
        return

    def trilegal_file_fmt(self):
        """ file name formats for trilegal input and trilegal sfr-z """
        tmp = os.path.split(self.cmd_input_file)[1]
        agb_mod = tmp.lower().replace('.dat', '').replace('cmd_input_', '')
        # trilegal output format
        self.tname = \
            os.path.join(self.outfile_loc, 'output_%s_%s_%s_%s' % (self.target,
                                                                   self.filter1,
                                                                   self.filter2,
                                                                   agb_mod))
        self.triout_fmt = self.tname + '_%003i.dat'

        sfr_fmt = '{}{}_tri_%003i.sfr'.format(self.target, self.filter1)
        self.sfr_fmt = os.path.join(self.outfile_loc, sfr_fmt)
        
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

    def prepare_trilegal_files(self, random_sfr=True, random_z=False,
                               zdisp=False, overwrite=False, object_mass=None):
        '''make the sfhs, make the galaxy inputs'''
        self.sfr_files = self.make_many_trilegal_sfhs(nsfhs=self.nsfhs,
                                                      outfile_fmt=self.sfr_fmt,
                                                      random_sfr=random_sfr,
                                                      random_z=random_z,
                                                      zdisp=zdisp,
                                                      dry_run=overwrite)

        self.prepare_galaxy_input(overwrite=overwrite, object_mass=object_mass)
        return

    def run_once(self, galaxy_input=None, triout=None, ite=0, overwrite=False):
        """write call to trilegal string"""
        flag = 0
        ver = 2.3
        call = ''

        #if os.path.isfile(triout) and not overwrite:
        #    logger.warning('{} exists, will overwrite if no hdf5 file found'.format(triout))
        #    flag += 1

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
        if self.nsfhs <= 1:
            cmd = self.run_once(galaxy_input=self.galaxy_input,
                                triout=self.tname + '_bestsfr.dat',
                                overwrite=overwrite)
        else:
            cmd = self.run_many(nproc=nproc, overwrite=overwrite)
        return cmd

    def run_many(self, nproc=8, overwrite=False):
        """Call self.run_once self.nsfh of times in iterations base on nproc"""
        self.prepare_trilegal_files(random_sfr=True, random_z=False,
                                    zdisp=False, overwrite=overwrite,
                                    object_mass=None)

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
                line += '{} &\n'.format(cmd)
            line += jobwait()
        return line

def call_VarySFH(input_file, loud=False, nproc=8, outfile=None,
                 overwrite=False):
    """
    write a script to run trilegal in parallel sampling the sfh
    
    overwrite may not work as expected:
    a) for trilegal input files (*.galinp, *.sfr): won't overwrite
    b) trilegal output: won't write the command to call trilegal if the .hdf5
       file already exists. (Will still overwrite the .dat file)
    """
    # set up logging
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

    vsh = VarySFHs(input_file=input_file)

    line = vsh.call_run(nproc=nproc, overwrite=overwrite)
    if outfile is None:
        print(line)
    else:
        logger.info('output file: {}'.format(outfile))
        with open(outfile, 'a') as out:
            out.write(line)
    return

def main(argv):
    """main function to call_VarySFH"""
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

    parser.add_argument('name', type=str,
                        help='input file. To create, see initialize_inputs.__doc__')

    args = parser.parse_args(argv)

    call_VarySFH(args.name, loud=args.verbose, nproc=args.nproc,
                 outfile=args.outfile, overwrite=args.overwrite)


if __name__ == '__main__':
    main(sys.argv[1:])
