import ResolvedStellarPops as rsp
import matplotlib.pyplot as plt
import fileIO
import matplotlib.transforms as mtransforms
try:
    from mpl_toolkits.axes_grid1.parasite_axes import SubplotHost
except:
    from mpl_toolkits.axes_grid.parasite_axes import SubplotHost
import os
from optparse import OptionParser
import itertools
import numpy as np
import shutil

import multiprocessing
import logging
logger = logging.getLogger()

sim_z = [0.0005, 0.001, 0.004, 0.008]

sim_age = [0.2, 1.0, 2.0, 8.0]

colorC = (0.8, 0.0, 0.0)
colorO = (0.0, 0.6, 0.8)


def ztomh(v):
    return np.log10(v / 0.019)


def write_sfh(age, z, ofile):
    oo = open(ofile, 'w')
    fmt = '%.6e %.2f %.3e\n'
    oo.write(fmt % (age - 1e-6, 0.0, z))
    oo.write(fmt % (age, 1.0, z))
    oo.write(fmt % (age + 1e-3, 1.0, z))
    oo.write(fmt % (age + 1e-3 + 1e-6, 0.0, z))
    oo.close()


def make_galaxy_input(galaxy_input, sfh_file):

    gal_dict_inp = {'photsys': '2mass',
                    'filter1': 'Ks',
                    'object_mass': 1e7,
                    'object_sfr_file': sfh_file}

    gal_dict = rsp.TrilegalUtils.galaxy_input_dict(**gal_dict_inp)

    gal_inp = rsp.fileIO.input_parameters(default_dict=gal_dict)

    galinp_kw = {'mag_limit_val': -3,
                 'object_sfr_mult_factorA': 1e9,
                 'mag_num': 3}

    gal_inp.add_params(galinp_kw)

    gal_inp.write_params(galaxy_input, rsp.TrilegalUtils.galaxy_input_fmt())
    return


def plot_em(ifile, lf_file, cmd_file, age, z, track):
    print 'Plotting', ifile

    logL, logTe, mbol, j, k, mcore, co, dmdt = np.loadtxt(ifile,
                                             usecols=(4, 5, 10, 11, 13, 14, 15, 18),
                                             unpack=True)

    nAGB = (mcore == 0)
    cAGB = ((co >= 1) & (dmdt <= -5))
    oAGB = ((co <= 1) & (logL >= 3.3) & (dmdt < -5))

    jk = j - k
    bins = np.arange(-10, 20, 0.1)

    ###### HRD
    fig = plt.figure()
    ax = fig.add_axes([.1, .1, .8, .8])
    ax.plot(logTe[nAGB], logL[nAGB], '.k')
    ax.plot(logTe[cAGB], logL[cAGB], 'o', mfc='None', ms=5, mew=1, mec=colorC,
            alpha=0.3)
    ax.plot(logTe[oAGB], logL[oAGB], 'o', mfc='None', ms=5, mew=1, mec=colorO,
            alpha=0.3)

    ax.annotate('Age=%.2e' % age, (.7, .1), va='center',
                xycoords='axes fraction')
    ax.annotate('Z=%.2e' % z, (.7, .15), va='center',
                xycoords='axes fraction')
    ax.annotate('[M/H]=%.2f' % ztomh(z), (.7, .2), va='center',
                xycoords='axes fraction')

    ax.annotate(r'$%s$' % track.replace('_', '\ '), (.1, .9), va='center',
                xycoords='axes fraction')

    ax.set_xlim(ax.get_xlim()[::-1])
    #ax.set_ylim(-3.1, -9.5)
    ax.set_xlabel(r'$\log\ T_{\\eff}$')
    ax.set_ylabel(r'$\log L$')
    plt.savefig(cmd_file.replace('cmd', 'hrd'))


    ###### CMD
    fig = plt.figure()
    ax = fig.add_axes([.1, .1, .8, .8])
    ax.plot(jk[nAGB], k[nAGB], '.k')
    ax.plot(jk[cAGB], k[cAGB], 'o', mfc='None', ms=5, mew=1, mec=colorC,
            alpha=0.3)
    ax.plot(jk[oAGB], k[oAGB], 'o', mfc='None', ms=5, mew=1, mec=colorO,
            alpha=0.3)

    ax.annotate('Age=%.2e' % age, (.7, .1), va='center',
                xycoords='axes fraction')
    ax.annotate('Z=%.2e' % z, (.7, .15), va='center',
                xycoords='axes fraction')
    ax.annotate('[M/H]=%.2f' % ztomh(z), (.7, .2), va='center',
                xycoords='axes fraction')

    ax.annotate(r'$%s$' % track.replace('_', '\ '), (.1, .9), va='center',
                xycoords='axes fraction')

    ax.set_xlim(.1, 2.4)
    ax.set_ylim(-3.1, -9.5)
    ax.set_xlabel(r'$J-K$')
    ax.set_ylabel(r'$K$')
    plt.savefig(cmd_file)
    plt.close()
    ###### LF
    fig = plt.figure()

    ax1 = SubplotHost(fig, 2, 1, 1)
    ax2 = SubplotHost(fig, 2, 1, 2)
    fig.add_subplot(ax1)
    fig.add_subplot(ax2)
    aux_trans = mtransforms.Affine2D().scale(-2.5, 1.).translate(4.77, 0)

    ax_logl = ax1.twin(aux_trans)
    ax_logl.set_viewlim_mode('transform')

    if sum(cAGB):
        pdf, bins, patches = ax1.hist(mbol[cAGB], bins, histtype='stepfilled',
                                      color=colorC)
        ax1.set_ylim(0.01, pdf.max()*1.1)

    if sum(oAGB):
        pdf, bins, patches = ax2.hist(mbol[oAGB], bins, histtype='stepfilled',
                                      color=colorO)
        ax2.set_ylim(0.01, pdf.max() * 1.1)

    ax1.annotate('Age=%.2e' % age, (.7, .1), va='center',
                 xycoords='axes fraction')
    ax1.annotate('Z=%.2e' % z, (.7, .2), va='center',
                 xycoords='axes fraction')
    ax1.annotate('[M/H]=%.2f' % ztomh(z), (.7, .3), va='center',
                 xycoords='axes fraction')

    ax1.annotate('C-rich', (.7, .8), va='center',
                 xycoords='axes fraction')
    ax2.annotate('O-rich', (.7, .8), va='center',
                 xycoords='axes fraction')

    ax2.annotate(r'$%s$' % track.replace('_', '\ '), (.1, .8), va='center',
                 xycoords='axes fraction')

    for ax in [ax1, ax2]:
        ax.set_xlim(0.2, -7.2)
        ax.axis['left'].set_label('N')

    #ax1.set_xlabel('Log L/L_sun')
    ax1.set_ylabel(r'$N$')
    ax2.set_ylabel(r'$N$')

    ax_logl.axis['right'].major_ticklabels.set_visible(False)
    ax1.axis['bottom'].major_ticklabels.set_visible(False)
    ax2.axis['bottom'].set_label('M$_\mathsf{bol}$')

    ax_logl.annotate(r'$\log L/L_\mathsf{sun}$', (.5, 1.1),
                     xycoords='axes fraction')
    #ax_logl.axis['top'].set_label('log L/L$_{\odot}$')
    fig.subplots_adjust(left=.1, bottom=None, right=0.98, top=None,
                        wspace=0, hspace=0)
    plt.savefig(lf_file)
    plt.close()

def run_all(age, z, track_set, sfh_dir, tri_dir, plt_dir, over_write=False):
    sfh_file = '%s/sfh_Z%.2e_A%.2e.dat' % (sfh_dir, z, age)
    inp_file = '%s/trilegal_input_Z%.2e_A%.2e.dat' % (tri_dir, z, age)
    out_file = '%s/trilegal_output_Z%.2e_A%.2e.dat' % (tri_dir, z, age)

    lf_file = '%s/lf_Z%.2e_A%.2e.png' % (plt_dir, z, age)
    cmd_file = '%s/cmd_Z%.2e_A%.2e.png' % (plt_dir, z, age)
    #print lf_file, cmd_file
    f_lf = open('%s/list_lf.dat' % plt_dir, 'a')
    f_cmd = open('%s/list_cmd.dat' % plt_dir, 'a')
    f_lf.write(lf_file+'\n')
    f_cmd.write(cmd_file+'\n')
    f_lf.close()
    f_cmd.close()

    write_sfh(age, z, sfh_file)
    make_galaxy_input(inp_file, sfh_file)
    track_file = '../cmd_inputfiles/cmd_input_%s.dat' % track_set
    rsp.TrilegalUtils.run_trilegal(track_file, inp_file, out_file,
                                   rmfiles=False)
    plot_em(out_file, lf_file, cmd_file, age, z, track_set)


def main(track_set, sfh_dir, tri_dir, plt_dir, over_write=False,
         multi=True):

    if os.path.isdir(plt_dir):
        shutil.rmtree(plt_dir)
        print 'rm\'d directories', plt_dir

    for d in [sfh_dir, tri_dir, plt_dir]:
        if not os.path.isdir(d):
            os.makedirs(d)
            print 'made directories', d
    if multi is False:        
        for age, z in itertools.product(sim_age, sim_z):
            run_all(age, z, track_set, sfh_dir, tri_dir, plt_dir)
            plt.close('all')
    else:
        pool = multiprocessing.Pool()
        res = []

        for age, z in itertools.product(sim_age, sim_z):    
            res.append(pool.apply_async(run_all,
                                        (age, z, track_set, sfh_dir, tri_dir, plt_dir),
                                        ))

        for r in res:
            r.get()

    return


if __name__ == '__main__':
    parser = OptionParser()

    parser.usage = '%prog track_set [options]'
    parser.add_option('--nr', action='store_true', default=False,
                      help='do not remove trilegal files')

    parser.add_option('-i', action='store_true', default=False,
                      help='supplying AGBTrackUtils inputfile')

    (options, args) = parser.parse_args()

    if len(args) != 1:
        parser.error('wrong args number')

    if options.i:
        input_file = args[0]
        infile = fileIO.input_file(input_file)
        agb_mix = infile.agb_mix
        set_name = infile.set_name
        track_set = '_'.join((agb_mix, set_name))
        diagnostic_dir = infile.diagnostic_dir0
        tri_dir = os.path.join(diagnostic_dir, 'trilegal_files')
        sfh_dir = os.path.join(tri_dir, 'sfh')
        plt_dir = os.path.join(diagnostic_dir, agb_mix, set_name, track_set)
        main(track_set, sfh_dir, tri_dir, plt_dir,
             over_write=infile.over_write, multi=False)
    else:
        cwd = os.getcwd()
        track_set = args[0]
        sfh_dir = os.path.join(cwd, 'SFH')
        tri_dir = os.path.join(cwd, 'TRILEGAL_FILES')
        plt_dir = os.path.join(cwd, 'PLOTS/%s' % track_set)
        main(track_set, sfh_dir, tri_dir, plt_dir)

    if not options.nr:
        shutil.rmtree(sfh_dir)
        shutil.rmtree(tri_dir)

    print
    print 'DONE'
