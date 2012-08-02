import os
import socket
hostname = socket.gethostname()

if hostname.endswith('astro.washington.edu'):
    print 'better edit TPAGBparams.py...'    
else:
    tpcalib_dir = '/Users/phil/research/TP-AGBcalib/'
    
table_src = os.path.join(tpcalib_dir,'SNAP/tables/')
data_src = os.path.join(tpcalib_dir,'SNAP/data')
output_src = os.path.join(tpcalib_dir,'SNAP/AST')
plt_dir = os.path.join(tpcalib_dir,'SNAP/plots')
model_src = os.path.join(tpcalib_dir,'SNAP/models')
fits_src = os.path.join(data_src,'galaxies')