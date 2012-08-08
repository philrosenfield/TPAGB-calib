import os
import socket
hostname = socket.gethostname()

if hostname.endswith('astro.washington.edu'):
    print 'better edit TPAGBparams.py...'    
else:
    tpcalib_dir = '/Users/phil/research/TP-AGBcalib/'

snap_src = os.path.join(tpcalib_dir,'SNAP')

table_src = os.path.join(snap_src,'tables')
plt_dir = os.path.join(snap_src,'plots','0.3dex/')
model_src = os.path.join(snap_src,'models','0.3dex/')
data_src = os.path.join(snap_src,'data')

output_src = os.path.join(model_src,'ast')
fits_src = os.path.join(data_src,'galaxies')