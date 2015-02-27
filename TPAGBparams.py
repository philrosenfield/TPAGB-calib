import os
import socket
hostname = socket.gethostname()

if hostname.endswith('astro.washington.edu'):
    print 'better edit TPAGBparams.py...'
else:
    if 'Linux' in os.uname():
        # linux laptop
        research_path = '/home/phil/research/'
        if os.uname()[1] == 'andromeda':
            # unipd
            research_path = '/home/rosenfield/research/'
    else:
        # mac
        research_path = '/Users/rosenfield/research/'

tpcalib_dir = research_path + 'TP-AGBcalib'

snap_src = os.path.join(tpcalib_dir, 'SNAP')
phat_src = os.path.join(tpcalib_dir, 'PHAT')

for directory in [snap_src, phat_src]:
    if not os.path.isdir(directory):
        print 'Warning. Can not find %s.' % directory
