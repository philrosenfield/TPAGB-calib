import sfh_tests
import os
import numpy as np
import matplotlib.pylab as plt


def load_plot_limits(filename='default'):
    if filename == 'default':
        filename = '/home/phil/research/TP-AGBcalib/SNAP/tables/cmd_plot_limits.dat'
    dtype = [('target', '|S16'),
             ('opt_cmdmin', '<f8'),
             ('opt_cmdmax', '<f8'),
             ('opt_lfmin', '<f8'),
             ('opt_lfmax', '<f8'),
             ('ir_cmdmin', '<f8'),
             ('ir_cmdmax', '<f8'),
             ('ir_lfmin', '<f8'),
             ('ir_lfmax', '<f8'),
             ('opt_offset', '<f8'),
             ('ir_offset', '<f8')]
    lims = np.genfromtxt(filename, dtype=dtype)
    return lims.view(np.recarray)


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
    lims = load_plot_limits()
    ax1.set_xlim(lims['target'==target]['opt_xmin'], lims['target'==target]['opt_xmax'])
    ax2.set_xlim(lims['target'==target]['ir_xmin'], lims['target'==target]['ir_xmax'])
    ax1.set_ylim(lims['target'==target]['opt_ymin'], lims['target'==target]['opt_ymax'])
    ax2.set_ylim(lims['target'==target]['ir_ymin'], lims['target'==target]['ir_ymax'])
    
    return ax1, ax2


def plot_lfs():
    pl = sfh_tests.Plotting()
    outfile_loc='/home/phil/Dropbox/research/varysfh/'
    cmd_inputs = ['CAF09_S_NOV13'.lower(),
                  'CAF09_S_NOV13eta0'.lower(),
                  'CAF09_S_OCT13'.lower()]
    targets = ['ddo71', 'hs117', 'kkh37', 'ngc2976-deep', 'ngc404', 'ddo78']
    one_plot = True
    lims = load_plot_limits()
    for target in targets:
        print target
        if one_plot is True:
            fig, axs = plt.subplots(ncols=2, figsize=(12,6))
            plt.subplots_adjust(right=0.95, left=0.05, wspace=0.1)
            cols = ['black', 'navy', 'darkgreen']
            narratio=False
        else:
            cols = ['black'] * len(cmd_inputs)
            axs = None
            narratio=True
        plot_data = True
        for i, cmd_input in enumerate(cmd_inputs):
            if i > 0 and one_plot is True:
                plot_data = False
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
                                         extra_str=cmd_input.split('_')[-1]+'_',
                                         axs=axs, plt_kw={'color': cols[i]},
                                         narratio=narratio, plot_data=plot_data)
            #ax1.set_title(cmd_input.replace('_', '\ '))
            lab = cmd_input.split('_')[-1]
            [ax.plot([0,0], [0,0], lw=3, color=cols[i], label='$%s$' % lab) for ax in [ax1, ax2]]
            row = lims[lims['target'] == target]
            
            ax1.set_xlim(row['opt_xmin'], row['opt_xmax'])
            ax2.set_xlim(row['ir_xmin'], row['ir_xmax'])
            ax1.set_ylim(row['opt_ymin'], row['opt_ymax'])
            ax2.set_ylim(row['ir_ymin'], row['ir_ymax'])
            [ax.legend(loc=0, frameon=False) for ax in [ax1, ax2]]
            figtitle = '%s%s_lfs.png' % (cmd_input.split('_')[-1]+'_', target)
            outfile = os.path.join(outfile_loc, figtitle)
