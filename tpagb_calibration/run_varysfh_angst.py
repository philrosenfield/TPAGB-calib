import os

from sfhs.vary_sfh import call_VarySFH
from TPAGBparams import snap_src

extra_strs = ['', '_555', '_475', '_606']

galaxy_lists = [['eso540-030',
                 'scl-de1',
                 'ugc-04459',
                 'ugc-4305-2',
                 'ugc-5139',
                 'ugc8508',
                 'ddo78',
                 'kdg73',
                 'kkh37',
                 'ngc3741',
                 'hs117',
                 'ngc2403-deep'],
                ['ic2574-sgs'],
                ['ugca292', 'ngc300-wide1'],
                ['ugca292', 'ngc300-wide1']]

vsfhpath = os.path.join(snap_src, 'varysfh')
for extra_str, galaxies in zip(extra_strs, galaxy_lists):
    for galaxy in galaxies:
        print galaxy
        input_file = os.path.join(vsfhpath, galaxy,
                                  '{}{}.vsfhinp'.format(galaxy, extra_str))

        call_VarySFH(input_file, loud=True, dry_run=False, max_proc=8)
