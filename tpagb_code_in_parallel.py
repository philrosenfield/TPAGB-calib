from IPython import parallel
import sfh_tests_multi_proc
import time
import numpy as np

def caller(vSFH, vsfh_kws):
    return vSFH.vary_the_SFH(**vsfh_kws)

clients = parallel.Client()
clients.block = False

clients[:].execute('cd ~/research/TP-AGBcalib/code/TPAGB-calib/')
clients[:].execute('import sfh_tests_multi_proc')

#reload(sfh_tests_multi_proc) 
#clients[:].execute('reload(sfh_tests_multi_proc)')

#targets = ['scl-de1', 'ddo71', 'hs117', 'kkh37', 'ngc2976-deep', 'ddo78']
#targets = ['ngc404']
targets = ['hs117', 'ngc2976-deep']
#cmd_inputs = ['cmd_input_CAF09_S_NOV13.dat', 'cmd_input_CAF09_S_NOV13eta0.dat', 'cmd_input_CAF09_S_OCT13.dat']
cmd_inputs = ['cmd_input_CAF09_S_FEB14.dat']
nsfhs = 15

vSFHs, vsfh_kws = sfh_tests_multi_proc.prepare_vsfh_run(targets, cmd_inputs, nsfhs)

# find a better way to run this all at once, rather then need to do to steps when nvsfhs > nprocs.
nprocs = len(clients)
nvsfhs = len(vSFHs)
ntimes = np.min([nprocs, nvsfhs])
ndiff = np.abs(nvsfhs - nprocs)

res = [clients[i].apply(caller, vSFHs[i], vsfh_kws,) for i in range(ntimes)]
while False in [r.ready() for r in res]:
    time.sleep(900)
[vSFHs[i].write_results(res[i].result) for i in range(ntimes)]
res2 = [clients[i].apply(caller, vSFHs[i+ntimes], vsfh_kws,) for i in range(ndiff)] 

while False in [r.ready() for r in res2]:
    time.sleep(900)
[vSFHs[i+ntimes].write_results(res2[i].result) for i in range(ndiff)]



