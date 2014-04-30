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
        research_path = '/Users/phil/research/'

tpcalib_dir = research_path + 'TP-AGBcalib'

snap_src = os.path.join(tpcalib_dir, 'SNAP')
phat_src = os.path.join(tpcalib_dir, 'PHAT')


table_src = os.path.join(snap_src, 'tables')
plt_dir = os.path.join(snap_src, 'plots')
model_src = os.path.join(snap_src, 'models')
data_src = os.path.join(snap_src, 'data')

output_src = os.path.join(model_src, 'ast')
fits_src = os.path.join(data_src, 'galaxies')

for directory in [research_path, tpcalib_dir, snap_src, table_src, plt_dir,
                  model_src, model_src, data_src, output_src, fits_src]:
    if not os.path.isdir(directory):
        print 'Warning. Can not find %s.' % directory
