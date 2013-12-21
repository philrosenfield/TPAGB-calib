import sfh_tests
import os

def plot_lf_with_stages(target, trilegal_output):
    pl = sfh_tests.Plotting()
    cols = ['navy', 'darkgreen', 'darkgreen', 'darkgreen', 'purple',
            'darkred']
    kw = {'outfile_loc': '/home/phil/Dropbox/research/varysfh/',
          'trilegal_output': trilegal_output,
          'narratio': False,
          'hist_it_up': False,
          'add_stage_lfs': 'all',
          'plot_data': False,
          'cols': cols,
          'stage_lf_kw': {'lw': 3}}
    ax1, ax2 = pl.compare_to_gal(target, **kw)
    ax1.set_xlim(23, ax1.get_xlim()[1])
    ax2.set_xlim(20.5, ax2.get_xlim()[1])
    ax1.set_ylim(10, 20000)
    ax2.set_ylim(10, 2000)
    
    return ax1, ax2

def plot_lfs():
    pl = sfh_tests.Plotting()
    outfile_loc='/home/phil/Dropbox/research/varysfh/'
    cmd_inputs = ['CAF09_S_NOV13'.lower(),
                  'CAF09_S_NOV13eta0'.lower(),
                  'CAF09_S_OCT13'.lower()]
    targets = ['ddo78', 'ddo71', 'hs117', 'kkh37', 'ngc2976-deep', 'ngc404-deep']
    for target in targets:
        for cmd_input in cmd_inputs:
            narratio_file_name = os.path.join(outfile_loc,
                                              '%s_%s_narratio.dat' %
                                              (cmd_input, target.lower()))
            opt_lf_file = os.path.join(outfile_loc,
                                       '%s_%s_opt_lf.dat' %
                                       (cmd_input, target.lower()))
            ir_lf_file =  os.path.join(outfile_loc,
                                       '%s_%s_ir_lf.dat' %
                                       (cmd_input, target.lower()))
        
            ax1, ax2 = pl.compare_to_gal(target, opt_lf_file=opt_lf_file,
                                         ir_lf_file=ir_lf_file,
                                         hist_it_up=False, outfile_loc=outfile_loc,
                                         narratio_file_name=narratio_file_name,
                                         extra_str=cmd_input.split('_')[-1]+'_')
            ax1.set_title(cmd_input.replace('_', '\ '))
