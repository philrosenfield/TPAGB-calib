from IPython import parallel
import brewer2mpl
import matplotlib.pylab as plt
import numpy as np
import os
import ResolvedStellarPops as rsp
import ResolvedStellarPops.graphics.GraphicsUtils as rg
import sys
import time


def caller(z, cmd_input, photsys, nsims=1, lagemin=7.1, lagemax=10.15,
           dlage=0.15):
    import os
    zdir = os.path.join(os.environ['CMDROOT'], 'MC_Z%.4f' % z)
    rsp.fileIO.ensure_dir(zdir)
    print 'Working on %s..' % zdir
    for i in range(nsims):
        fname = 'caf09_s_nov13_%s_Z%.4f_%003d.dat' % (photsys, z, i)
        print '%i/%s: %s' % ((i+1), nsims, fname)
        outfile = os.path.join(zdir, fname)
        if os.path.isfile(outfile):
            continue
        rsp.IsochroneUtils.run_isoch(cmd_input, outfile, photsys,
               isoch_input_kw={'other_qty': z, 'qty_min': lagemin,
                               'qty_max': lagemax, 'dq': dlage,
                               'kind_iso': 5},
               cmd_input_kw={'file_imf': 'tab_imf/imf_kroupa.dat'})
    return


def integrated_colors(cmd_input, photsys, zs='default', dry_run=False):
    here = os.getcwd()
    clients = parallel.Client()
    clients.block = False

    clients[:].execute('cd $CMDROOT')
    clients[:].execute('import ResolvedStellarPops as rsp')
    clients[:].execute('import os')

    if dry_run is True:
        timeout = 10
    else:
        timeout = 30

    if zs == 'default':
        zs = [0.0005, 0.001, 0.002, 0.003, 0.004, 0.005, 0.006, 0.007, 0.008,
              0.01, 0.02]

    # find a better way to run this all at once, what if I need three times
    # through?
    nprocs = len(clients)
    nneeded = len(zs)
    ntimes = np.min([nprocs, nneeded])
    ndiff = np.abs(nneeded - nprocs)
    print nprocs, nneeded, ntimes, ndiff
    if ndiff > nprocs:
        print 'need a for loop, too many processes code code code man'
        import sys
        sys.exit()

    print 'calling first set'
    res = [clients[i].apply(caller, zs[i], cmd_input, photsys,)
           for i in range(ntimes)]
    #res = [caller(zs[i], cmd_input, photsys) for i in range(nneeded)]
    #sys.exit()
    while False in [r.ready() for r in res]:
        print 'waiting on first set'
        time.sleep(timeout)
        print 'checking first set...'

    if not nprocs > nneeded:
        print 'calling second set'
        res2 = [clients[i].apply(caller, zs[i+ntimes], cmd_input, photsys,)
                for i in range(ndiff)]

        while False in [r.ready() for r in res2]:
            print 'waiting on second set'
            time.sleep(timeout)
            print 'checking second set...'

    print 'done.'
    os.chdir(here)

def plot_integrated_colors(filenames, labels='Z'):
    if type(filenames) is str:
        filenames = [filenames]
        ax = None
        cols = ['k']
    else:
        fig, ax = plt.subplots()
        cols = brewer2mpl.get_map('Spectral', 'Diverging',
                                  len(filenames)).mpl_colors

    if labels == 'Z':
        fmt = '$Z=%.4f$'
        labels = [fmt % float(l.replace('.dat', '').split('Z')[1])
                  for l in filenames]
    else:
        print 'need to fix labels'
        labels = [''] * len(filenames)
    for i, filename in enumerate(filenames):
        data = rsp.fileIO.readfile(filename)
        ycol = 'V-K'
        xcol = 'Age'
        ax = rg.color_color(data, xcol, ycol, xscale='log', ax=ax,
                            plt_kw={'lw': 2, 'color': cols[i],
                                    'label': labels[i]})

    plot_cluster_data(ax)
    ax.legend(frameon=False, loc=0, numpoints=1)
    ax.set_xlabel(r'${\rm %s}$' % xcol, fontsize=20)
    ax.set_ylabel(r'${\rm %s}$' % ycol, fontsize=20)
    plt.tick_params(labelsize=16)
    return ax

def plot_mc_integrated_colors(filenames):
    fig, ax = plt.subplots()
    for i, filename in enumerate(filenames):
        data = rsp.fileIO.readfile(filename)
        ycol = 'V-K'
        xcol = 'Age'
        ax = rg.color_color(data, xcol, ycol, xscale='log', ax=ax,
                            plt_kw={'marker': '.', 'color': 'k',
                                    'alpha': 0.3})

    plot_cluster_data(ax)

def s_to_age(s):
    s[np.abs(s)==100] = 53.03
    return 10**(6.227 + 0.0733 * s)

def plot_cluster_data(ax=None):
    if ax is None:
        fig, ax = plt.subplots()
    lmc_age, lmc_vk, kyeong_age, kyeong_vk, pas_age, pas_vk = cluster_data()
    for age, vk, lab, mark, col in zip([lmc_age, kyeong_age, pas_age],
                       [lmc_vk, kyeong_vk, pas_vk],
                       ['P83', 'K03', 'G06 + P06'],
                       ['o', 's', '*'],
                       ['darkgreen', 'navy', 'darkred']):
        ax.scatter(age, vk, marker=mark, label=lab, color=col)
    return ax


def cluster_data():
    loc = '/home/rosenfield/research/TP-AGBcalib/clusters'
    ir_lmc = rsp.fileIO.readfile(os.path.join(loc, 'ir_lmc.dat'))
    inds, = np.nonzero(ir_lmc['VK'] > 0)
    lmc_age = s_to_age(ir_lmc['S'][inds])
    lmc_vk = ir_lmc['VK'][inds]

    kyeong = rsp.fileIO.readfile(os.path.join(loc, 'Kyeong_VK.dat'), col_key_line=2)
    kyeong_age = s_to_age(kyeong['S'])

    kyeong_vk = kyeong['V'] - kyeong['K']

    pas = rsp.fileIO.readfile(os.path.join(loc, 'Pas_VK_ok.dat'), col_key_line=1)
    inds, = np.nonzero(pas['V'] > 0)
    pas_age = s_to_age(pas['sL'][inds])
    pas_vk = pas['V'][inds] - pas['K'][inds]

    return lmc_age, lmc_vk, kyeong_age, kyeong_vk, pas_age, pas_vk



if __name__ == '__main__':

    args = sys.argv
    if len(args) == 1:
        cmd_input = 'cmd_input_parsecCAF09_S12_NS_1TP_NOV13.dat'
        photsys = 'ubvrijhk'
        nproc = 12
        print 'using default values'
    else:
        cmd_input, photsys, nproc = sys.argv[1:]
    print 'cmd_input:', cmd_input
    print 'photsys:', photsys
    os.system('ipcluster start -n=%i &' % nproc)
    timeout = nproc * 5
    print 'waiting %i sec for clusters to spin up...' % timeout
    time.sleep(timeout)
    print 'moving on'

    integrated_colors(cmd_input, photsys)
    os.system('ipcluster stop &')
